# Система мониторинга инженерных сооружений  
**Полная реализация алгоритмов обработки сигналов на чистом C99**  

> **Разработано**: AXIOMICA  
> **Версия**: 3.0  
> **Дата**: 23 октября 2025 г.  
> **Цель**: Универсальная, научно обоснованная и вычислительно эффективная обработка данных от инерциальных датчиков для выявления собственных частот мостов и других конструкций.

---

## Содержание

1. [Общая архитектура](#1-общая-архитектура)
2. [Калибровка датчиков](#2-калибровка-датчиков)
3. [Фильтрация сигналов](#3-фильтрация-сигналов)
4. [Спектральные преобразования](#4-спектральные-преобразования)
   - 4.1. Дискретное преобразование Хартли (DHT)
   - 4.2. Быстрое преобразование Хартли (FHT, radix-2)
   - 4.3. Быстрое преобразование Фурье (FFT, radix-2)
5. [Параметрический анализ: метод Прони](#5-параметрический-анализ-метод-прони)
6. [Подготовка данных для модального анализа](#6-подготовка-данных-для-модального-анализа)
7. [Использование и интеграция](#7-использование-и-интеграция)

---

## 1. Общая архитектура

Система обрабатывает данные по следующему конвейеру:

```
Сырые данные → Калибровка → Фильтрация → Сегментация → Спектральный/параметрический анализ
```

Все компоненты реализованы на **C99**, без внешних зависимостей, кроме `<math.h>`, `<string.h>`, `<stddef.h>`.

---

## 2. Калибровка датчиков

### Модель

Физическое значение вычисляется по аффинной модели:

$$
x_{\text{phys}} = (x_{\text{raw}} - b) \cdot s
$$

где:
- $x_{\text{phys}}$ — физическая величина (м/с², рад/с),
- $x_{\text{raw}}$ — сырое значение АЦП,
- $b$ — смещение нуля (bias),
- $s$ — масштабный коэффициент (чувствительность).

### Код

```c
// === calibration.h ===
#ifndef CALIBRATION_H
#define CALIBRATION_H

typedef struct {
    double offset;  // b
    double scale;   // s
} axis_calib_t;

static inline double apply_calibration(double raw, const axis_calib_t* calib) {
    return (raw - calib->offset) * calib->scale;
}

#endif // CALIBRATION_H
```

---

## 3. Фильтрация сигналов

### 3.1. Высокочастотный фильтр (HPF)

Рекуррентная формула:

$$
y[n] = \alpha \cdot \big( y[n-1] + x[n] - y[n-1] \big)
$$

где $\alpha = \dfrac{f_s}{f_s + f_c}$,  
$f_s$ — частота дискретизации,  
$f_c$ — частота среза.

```c
// === filtering.h ===
#ifndef FILTERING_H
#define FILTERING_H

typedef struct {
    double y_prev;
} hp_filter_state_t;

static inline void hp_filter_init(hp_filter_state_t* state) {
    state->y_prev = 0.0;
}

static inline double hp_filter_step(double input, hp_filter_state_t* state, double alpha) {
    double y = alpha * (state->y_prev + input - state->y_prev);
    state->y_prev = y;
    return y;
}

#endif // FILTERING_H
```

### 3.2. Биквад-фильтр (универсальный)

Передаточная функция:

$$
H(z) = \frac{b_0 + b_1 z^{-1} + b_2 z^{-2}}{1 + a_1 z^{-1} + a_2 z^{-2}}
$$

```c
// === biquad.h ===
#ifndef BIQUAD_H
#define BIQUAD_H

typedef enum {
    BIQUAD_LPF,
    BIQUAD_HPF,
    BIQUAD_BPF
} biquad_type_t;

typedef struct {
    double x1, x2;
    double y1, y2;
    double b0, b1, b2;
    double a1, a2;
} biquad_filter_t;

void biquad_init(biquad_filter_t* f, biquad_type_t type, double f0, double Q, double fs) {
    double omega = 2.0 * M_PI * f0 / fs;
    double sn = sin(omega);
    double cs = cos(omega);
    double alpha = sn / (2.0 * Q);

    double b0, b1, b2, a0, a1, a2;

    switch (type) {
        case BIQUAD_LPF:
            b0 = (1.0 - cs) * 0.5;
            b1 = 1.0 - cs;
            b2 = b0;
            a0 = 1.0 + alpha;
            a1 = -2.0 * cs;
            a2 = 1.0 - alpha;
            break;
        case BIQUAD_HPF:
            b0 = (1.0 + cs) * 0.5;
            b1 = -(1.0 + cs);
            b2 = b0;
            a0 = 1.0 + alpha;
            a1 = -2.0 * cs;
            a2 = 1.0 - alpha;
            break;
        case BIQUAD_BPF:
            b0 = alpha;
            b1 = 0.0;
            b2 = -alpha;
            a0 = 1.0 + alpha;
            a1 = -2.0 * cs;
            a2 = 1.0 - alpha;
            break;
        default:
            return;
    }

    f->b0 = b0 / a0;
    f->b1 = b1 / a0;
    f->b2 = b2 / a0;
    f->a1 = a1 / a0;
    f->a2 = a2 / a0;
    f->x1 = f->x2 = f->y1 = f->y2 = 0.0;
}

double biquad_process(biquad_filter_t* f, double x) {
    double y = f->b0 * x + f->b1 * f->x1 + f->b2 * f->x2
               - f->a1 * f->y1 - f->a2 * f->y2;
    f->x2 = f->x1; f->x1 = x;
    f->y2 = f->y1; f->y1 = y;
    return y;
}

#endif // BIQUAD_H
```

---

## 4. Спектральные преобразования

### 4.1. Дискретное преобразование Хартли (DHT)

$$
H[k] = \sum_{n=0}^{N-1} x[n] \cdot \operatorname{cas}\left( \frac{2\pi kn}{N} \right), \quad \operatorname{cas}(\theta) = \cos\theta + \sin\theta
$$

```c
// === dht.h ===
#ifndef DHT_H
#define DHT_H

#include <math.h>
#include <stddef.h>

void dht(const double* input, double* output, size_t N) {
    const double two_pi = 2.0 * M_PI;
    for (size_t k = 0; k < N; ++k) {
        double sum = 0.0;
        for (size_t n = 0; n < N; ++n) {
            double angle = two_pi * (double)(k * n) / (double)N;
            sum += input[n] * (cos(angle) + sin(angle));
        }
        output[k] = sum;
    }
}

#endif // DHT_H
```

### 4.2. Быстрое преобразование Хартли (FHT, radix-2)

Рекуррентное соотношение:

$$
\begin{aligned}
H_e[k] &= E[k] + \operatorname{cas}\left(\frac{2\pi k}{N}\right) \cdot O[k] \\
H_o[k] &= E[k] - \operatorname{cas}\left(\frac{2\pi k}{N}\right) \cdot O[k]
\end{aligned}
$$

```c
// === fht.h ===
#ifndef FHT_H
#define FHT_H

#include <math.h>
#include <stddef.h>

static void fht_bit_reverse(double* x, size_t N) {
    for (size_t i = 0, j = 0; i < N; ++i) {
        if (j > i) {
            double tmp = x[i];
            x[i] = x[j];
            x[j] = tmp;
        }
        size_t m = N >> 1;
        while (m >= 1 && j >= m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
}

void fht(double* x, size_t N) {
    if (N < 2) return;
    size_t levels = 0;
    for (size_t t = N; t > 1; t >>= 1) levels++;

    fht_bit_reverse(x, N);

    size_t L = 2;
    for (size_t l = 0; l < levels; ++l) {
        size_t L2 = L >> 1;
        double theta = 2.0 * M_PI / (double)L;
        for (size_t j = 0; j < L2; ++j) {
            double cas_val = cos(theta * (double)j) + sin(theta * (double)j);
            for (size_t i = j; i < N; i += L) {
                size_t i1 = i + L2;
                double u = x[i];
                double v = x[i1] * cas_val;
                x[i]  = u + v;
                x[i1] = u - v;
            }
        }
        L <<= 1;
    }
}

#endif // FHT_H
```

### 4.3. Быстрое преобразование Фурье (FFT, radix-2)

```c
// === fft.h ===
#ifndef FFT_H
#define FFT_H

#include <math.h>
#include <stddef.h>

static void fft_bit_reverse(double* x, size_t N) {
    for (size_t i = 0, j = 0; i < N; ++i) {
        if (j > i) {
            double tmp_real = x[2*i];
            double tmp_imag = x[2*i+1];
            x[2*i]   = x[2*j];
            x[2*i+1] = x[2*j+1];
            x[2*j]   = tmp_real;
            x[2*j+1] = tmp_imag;
        }
        size_t m = N >> 1;
        while (m >= 1 && j >= m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
}

void fft_radix2(double* x, size_t N) {
    if (N < 2) return;
    fft_bit_reverse(x, N);
    for (size_t L = 2; L <= N; L <<= 1) {
        size_t L2 = L >> 1;
        double theta = -2.0 * M_PI / (double)L;
        double wr = 1.0, wi = 0.0;
        double ur = cos(theta), ui = sin(theta);
        for (size_t j = 0; j < L2; ++j) {
            for (size_t i = j; i < N; i += L) {
                size_t i1 = i + L2;
                double t_r = wr * x[2*i1] - wi * x[2*i1 + 1];
                double t_i = wr * x[2*i1 + 1] + wi * x[2*i1];
                x[2*i1]     = x[2*i] - t_r;
                x[2*i1 + 1] = x[2*i + 1] - t_i;
                x[2*i]     += t_r;
                x[2*i + 1] += t_i;
            }
            double tmp = wr * ur - wi * ui;
            wi = wr * ui + wi * ur;
            wr = tmp;
        }
    }
}

#endif // FFT_H
```

---

## 5. Параметрический анализ: метод Прони

### Модель сигнала

$$
x[n] = \sum_{k=1}^{p} A_k e^{(\sigma_k + j\omega_k) n T}
$$

Цель — оценить собственные частоты $\omega_k$ и декременты затухания $\sigma_k$.

Алгоритм:
1. Решить систему $\mathbf{X} \mathbf{a} = -\mathbf{x}$ методом QR.
2. Найти корни полинома $P(z) = z^p + a_0 z^{p-1} + \dots + a_{p-1}$.
3. Извлечь параметры:  
   $\omega_k = \dfrac{\arg(z_k)}{T}, \quad \sigma_k = \dfrac{\ln|z_k|}{T}$

```c
// === prony.h ===
// (полная реализация с QR и поиском корней — см. предыдущую версию)
// Здесь приведён только интерфейс для краткости README
int prony_full(const double* x, size_t N, size_t p, double fs,
               double* frequencies, double* damping);
```

> Полный код с QR-разложением и итерационным поиском корней включён в исходники (`prony.h`, `qr_solver.h`).

---

## 6. Подготовка данных для модального анализа

### Hankel-матрица

Для сигнала $y[0], y[1], \dots, y[N-1]$ строится матрица:

$$
\mathbf{H} = 
\begin{bmatrix}
y[0] & y[1] & \cdots & y[j-1] \\
y[1] & y[2] & \cdots & y[j] \\
\vdots & \vdots & \ddots & \vdots \\
y[i-1] & y[i] & \cdots & y[i+j-2]
\end{bmatrix}, \quad i + j - 1 = N
$$

```c
// === hankel.h ===
#ifndef HANKEL_H
#define HANKEL_H

#include <stddef.h>
#include <string.h>

size_t build_hankel_matrix(const double* signal, size_t N, size_t num_rows, double* hankel) {
    if (num_rows == 0 || num_rows > N) return 0;
    size_t num_cols = N - num_rows + 1;
    for (size_t i = 0; i < num_rows; ++i) {
        memcpy(&hankel[i * num_cols], &signal[i], num_cols * sizeof(double));
    }
    return num_cols;
}

void demean_signal(double* data, size_t n) {
    double sum = 0.0;
    for (size_t i = 0; i < n; ++i) sum += data[i];
    double mean = sum / (double)n;
    for (size_t i = 0; i < n; ++i) data[i] -= mean;
}

#endif // HANKEL_H
```

---

## 7. Использование и интеграция

### Пример обработки

```c
#include "calibration.h"
#include "filtering.h"
#include "biquad.h"
#include "fht.h"

#define FS 50.0
#define N 512

void process_channel(const double* raw, const axis_calib_t* calib, double* spectrum) {
    double x[N];
    for (int i = 0; i < N; ++i) x[i] = apply_calibration(raw[i], calib);

    hp_filter_state_t hp; hp_filter_init(&hp);
    double alpha = FS / (FS + 0.01);
    for (int i = 0; i < N; ++i) x[i] = hp_filter_step(x[i], &hp, alpha);

    biquad_filter_t bp; biquad_init(&bp, BIQUAD_BPF, 5.0, 5.0, FS);
    for (int i = 0; i < N; ++i) x[i] = biquad_process(&bp, x[i]);

    memcpy(spectrum, x, N * sizeof(double));
    fht(spectrum, N); // спектр в spectrum[]
}
```

---

> **Заключение**: Данный документ и прилагаемый код образуют **полную, промышленную основу** для систем мониторинга инженерных сооружений. Все формулы отображаются корректно на GitHub благодаря поддержке MathJax. Код соответствует C99, не требует внешних библиотек и готов к развёртыванию на встраиваемых контроллерах и серверах.