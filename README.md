```markdown
# Система мониторинга инженерных сооружений  
**Полная реализация алгоритмов обработки сигналов на чистом C99**

> **Разработано**: AXIOMICA  
> **Версия**: 4.0  
> **Дата**: 23 октября 2025 г.

---

## Содержание

1. [Общая архитектура](#1-общая-архитектура)
2. [Калибровка датчиков](#2-калибровка-датчиков)
3. [Фильтрация сигналов](#3-фильтрация-сигналов)
4. [Спектральные преобразования](#4-спектральные-преобразования)
5. [Параметрический анализ: метод Прони](#5-параметрический-анализ-метод-прони)
6. [Подготовка данных для модального анализа](#6-подготовка-данных-для-модального-анализа)
7. [Использование](#7-использование)

---

## 1. Общая архитектура

Конвейер обработки:

```
Сырые данные → Калибровка → Фильтрация → Сегментация → Анализ
```

Все компоненты — на C99, без внешних зависимостей, кроме `<math.h>`, `<string.h>`, `<stddef.h>`.

---

## 2. Калибровка датчиков

Модель:

$$
x_{\text{phys}} = (x_{\text{raw}} - b) \cdot s
$$

```c
// === calibration.h ===
#ifndef CALIBRATION_H
#define CALIBRATION_H

/**
 * @brief Параметры аффинной калибровки одной оси датчика.
 *
 * Используются для преобразования сырых значений АЦП в физические единицы.
 * Коэффициенты определяются в процессе лабораторной калибровки.
 */
typedef struct {
    double offset;  ///< Смещение нуля (bias), единицы АЦП.
    double scale;   ///< Масштабный коэффициент (чувствительность), физ.ед./АЦП.
} axis_calib_t;

/**
 * @brief Применяет аффинную калибровку к одному значению.
 *
 * @param raw   Сырое значение от датчика (АЦП).
 * @param calib Указатель на структуру с калибровочными коэффициентами.
 * @return      Значение в физических единицах.
 *
 * @note Функция inline для минимизации накладных расходов в циклах обработки.
 */
static inline double apply_calibration(double raw, const axis_calib_t* calib) {
    return (raw - calib->offset) * calib->scale;
}

#endif // CALIBRATION_H
```

---

## 3. Фильтрация сигналов

### 3.1. Высокочастотный фильтр

$$
y[n] = \alpha \cdot (y[n-1] + x[n] - y[n-1]), \quad \alpha = \frac{f_s}{f_s + f_c}
$$

```c
// === filtering.h ===
#ifndef FILTERING_H
#define FILTERING_H

/**
 * @brief Состояние рекурсивного высокочастотного фильтра первого порядка.
 *
 * Фильтр предназначен для подавления медленного дрейфа и постоянной составляющей.
 * Требует сохранения состояния между вызовами.
 */
typedef struct {
    double y_prev;  ///< Предыдущее выходное значение фильтра.
} hp_filter_state_t;

/**
 * @brief Инициализация состояния HPF.
 *
 * @param state Указатель на структуру состояния.
 *
 * @note Обнуляет предыдущее значение. Вызывать один раз перед началом фильтрации.
 */
static inline void hp_filter_init(hp_filter_state_t* state) {
    state->y_prev = 0.0;
}

/**
 * @brief Обработка одного отсчёта высокочастотным фильтром.
 *
 * @param input Входное значение (откалиброванное).
 * @param state Указатель на состояние фильтра (должен быть инициализирован).
 * @param alpha Коэффициент фильтра: alpha = fs / (fs + fc).
 * @return      Отфильтрованное значение.
 *
 * @note Коэффициент alpha должен быть предварительно рассчитан на основе fs и fc.
 */
static inline double hp_filter_step(double input, hp_filter_state_t* state, double alpha) {
    double y = alpha * (state->y_prev + input - state->y_prev);
    state->y_prev = y;
    return y;
}

#endif // FILTERING_H
```

### 3.2. Биквад-фильтр

$$
H(z) = \frac{b_0 + b_1 z^{-1} + b_2 z^{-2}}{1 + a_1 z^{-1} + a_2 z^{-2}}
$$

```c
// === biquad.h ===
#ifndef BIQUAD_H
#define BIQUAD_H

/**
 * @brief Типы биквад-фильтров.
 */
typedef enum {
    BIQUAD_LPF,  ///< Фильтр низких частот (Butterworth, 2-го порядка)
    BIQUAD_HPF,  ///< Фильтр высоких частот (Butterworth, 2-го порядка)
    BIQUAD_BPF   ///< Полосовой фильтр (Butterworth, 2-го порядка)
} biquad_type_t;

/**
 * @brief Состояние и коэффициенты биквад-фильтра.
 *
 * Реализует рекурсивный IIR-фильтр второго порядка в прямой форме II.
 */
typedef struct {
    double x1, x2;  ///< Задержанные входные значения (x[n-1], x[n-2])
    double y1, y2;  ///< Задержанные выходные значения (y[n-1], y[n-2])
    double b0, b1, b2;  ///< Коэффициенты числителя передаточной функции
    double a1, a2;      ///< Коэффициенты знаменателя (a0 = 1 подразумевается)
} biquad_filter_t;

/**
 * @brief Инициализация биквад-фильтра.
 *
 * @param f    Указатель на структуру фильтра.
 * @param type Тип фильтра (LPF, HPF, BPF).
 * @param f0   Центральная частота (для BPF) или частота среза (для LPF/HPF), Гц.
 * @param Q    Добротность. Для BPF: Q = f0 / BW. Для Butterworth: Q = 1/sqrt(2).
 * @param fs   Частота дискретизации, Гц.
 *
 * @note Используется билинейное преобразование с предварительной коррекцией частоты.
 */
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

/**
 * @brief Обработка одного отсчёта биквад-фильтром.
 *
 * @param f Указатель на инициализированную структуру фильтра.
 * @param x Входное значение.
 * @return  Отфильтрованное значение.
 *
 * @note Функция потокобезопасна только если структура f не разделяется между потоками.
 */
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

/**
 * @brief Прямое дискретное преобразование Хартли (DHT).
 *
 * @param input  Указатель на входной массив длины N (вещественные числа).
 * @param output Указатель на выходной массив длины N.
 * @param N      Длина последовательности (должна быть > 0).
 *
 * @note Алгоритм имеет сложность O(N^2). Использовать только для N <= 512.
 * @note Обратное преобразование совпадает с прямым с точностью до масштаба 1/N.
 */
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

```c
// === fht.h ===
#ifndef FHT_H
#define FHT_H

#include <math.h>
#include <stddef.h>

/**
 * @brief Бит-реверсивная перестановка (in-place).
 *
 * @param x Массив длины N (должен быть степенью двойки).
 * @param N Длина массива.
 */
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

/**
 * @brief Быстрое преобразование Хартли (FHT, radix-2, in-place).
 *
 * @param x Указатель на массив данных длины N (должен быть степенью двойки).
 * @param N Длина массива (степень двойки, N >= 2).
 *
 * @note Алгоритм основан на рекуррентном соотношении с cas-функцией.
 * @note Входной массив заменяется результатом преобразования.
 */
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

/**
 * @brief Бит-реверсивная перестановка для комплексного FFT.
 *
 * @param x Массив длины 2*N: [Re0, Im0, Re1, Im1, ..., Re{N-1}, Im{N-1}]
 * @param N Длина комплексной последовательности (степень двойки).
 */
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

/**
 * @brief Radix-2 FFT (in-place, decimation-in-time).
 *
 * @param x Указатель на массив длины 2*N: [Re0, Im0, ..., Re{N-1}, Im{N-1}]
 * @param N Длина комплексной последовательности (степень двойки, N >= 2).
 *
 * @note Алгоритм использует стандартную butterfly-операцию.
 * @note Для вещественного входа можно использовать оптимизацию real-FFT.
 */
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

Модель сигнала:

$$
x[n] = \sum_{k=1}^{p} A_k e^{(\sigma_k + j\omega_k) n T}
$$

```c
// === qr_solver.h ===
#ifndef QR_SOLVER_H
#define QR_SOLVER_H

#include <math.h>
#include <stddef.h>

#define MAX_PRONY_ORDER 16
#define MAX_PRONY_SAMPLES 1024

/**
 * @brief Решение системы линейных уравнений A * x = b методом QR-разложения (Householder).
 *
 * @param A Матрица коэффициентов размера m x n (m >= n), хранится по строкам.
 * @param b Вектор правой части длины m.
 * @param x Выходной вектор решения длины n.
 * @param m Число уравнений.
 * @param n Число неизвестных.
 * @return  0 при успехе, -1 при ошибке (вырожденная матрица или превышение размеров).
 *
 * @note Использует статические буферы для избежания динамического выделения памяти.
 * @note Подходит для задач небольшой размерности (n <= 16).
 */
int qr_solve(double* A, double* b, double* x, size_t m, size_t n) {
    if (n > MAX_PRONY_ORDER || m > MAX_PRONY_SAMPLES) return -1;
    if (m < n) return -1;

    double R[MAX_PRONY_ORDER * MAX_PRONY_ORDER];
    double Qtb[MAX_PRONY_ORDER];
    memcpy(R, A, m * n * sizeof(double));

    // Householder reflections
    for (size_t k = 0; k < n; ++k) {
        double norm = 0.0;
        for (size_t i = k; i < m; ++i) {
            norm += R[i * n + k] * R[i * n + k];
        }
        if (norm < 1e-15) return -1;
        double alpha = sqrt(norm);
        if (R[k * n + k] < 0) alpha = -alpha;
        R[k * n + k] += alpha;
        norm = 0.0;
        for (size_t i = k; i < m; ++i) {
            norm += R[i * n + k] * R[i * n + k];
        }
        norm = sqrt(norm);
        if (norm < 1e-15) continue;
        for (size_t i = k; i < m; ++i) {
            R[i * n + k] /= norm;
        }

        // Применяем отражение к A и b
        for (size_t j = k + 1; j < n; ++j) {
            double dot = 0.0;
            for (size_t i = k; i < m; ++i) {
                dot += R[i * n + k] * R[i * n + j];
            }
            for (size_t i = k; i < m; ++i) {
                R[i * n + j] -= 2.0 * dot * R[i * n + k];
            }
        }

        double dot_b = 0.0;
        for (size_t i = k; i < m; ++i) {
            dot_b += R[i * n + k] * b[i];
        }
        for (size_t i = k; i < m; ++i) {
            b[i] -= 2.0 * dot_b * R[i * n + k];
        }
    }

    // Решаем R(0:n,0:n) * x = b(0:n)
    for (size_t i = 0; i < n; ++i) {
        Qtb[i] = b[i];
    }

    for (ptrdiff_t i = (ptrdiff_t)n - 1; i >= 0; --i) {
        double sum = Qtb[i];
        for (size_t j = (size_t)(i + 1); j < n; ++j) {
            sum -= R[i * n + j] * x[j];
        }
        if (fabs(R[i * n + i]) < 1e-12) return -1;
        x[i] = sum / R[i * n + i];
    }
    return 0;
}

#endif // QR_SOLVER_H
```

```c
// === prony.h ===
#ifndef PRONY_H
#define PRONY_H

#include "qr_solver.h"
#include <math.h>
#include <string.h>

/**
 * @brief Полная реализация метода Прони для оценки модальных параметров.
 *
 * @param x           Входной сигнал длины N.
 * @param N           Длина сигнала.
 * @param p           Порядок модели (макс. число мод).
 * @param fs          Частота дискретизации, Гц.
 * @param frequencies Выходной массив частот (Гц), длина p.
 * @param damping     Выходной массив декрементов затухания, длина p.
 * @return            Число успешно оценённых мод (<= p).
 *
 * @note Алгоритм:
 *       1. Формирование системы уравнений Прони.
 *       2. Решение методом QR.
 *       3. Поиск корней характеристического полинома методом Ньютона.
 * @note Требует N >= 2*p.
 */
int prony_full(const double* x, size_t N, size_t p, double fs,
               double* frequencies, double* damping) {
    if (p == 0 || p > MAX_PRONY_ORDER || N < 2 * p) return 0;

    double A[MAX_PRONY_SAMPLES * MAX_PRONY_ORDER];
    double b[MAX_PRONY_SAMPLES];
    size_t rows = N - p;

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < p; ++j) {
            A[i * p + j] = x[i + p - 1 - j];
        }
        b[i] = -x[i + p];
    }

    double a[MAX_PRONY_ORDER];
    if (qr_solve(A, b, a, rows, p) != 0) return 0;

    int found = 0;
    for (size_t k = 0; k < p && found < (int)p; ++k) {
        double theta0 = 2.0 * M_PI * (double)k / (double)p;
        double z_re = cos(theta0);
        double z_im = sin(theta0);

        for (int iter = 0; iter < 50; ++iter) {
            double p_re = 1.0, p_im = 0.0;
            double dp_re = (double)p, dp_im = 0.0;

            for (size_t i = 0; i < p; ++i) {
                double new_re = p_re * z_re - p_im * z_im + a[i];
                double new_im = p_re * z_im + p_im * z_re;
                p_re = new_re;
                p_im = new_im;

                if (i < p - 1) {
                    double coef = (double)(p - 1 - i);
                    double new_dp_re = dp_re * z_re - dp_im * z_im + coef * a[i];
                    double new_dp_im = dp_re * z_im + dp_im * z_re;
                    dp_re = new_dp_re;
                    dp_im = new_dp_im;
                }
            }

            double denom = dp_re * dp_re + dp_im * dp_im;
            if (denom < 1e-15) break;

            double dz_re = (p_re * dp_re + p_im * dp_im) / denom;
            double dz_im = (p_im * dp_re - p_re * dp_im) / denom;

            z_re -= dz_re;
            z_im -= dz_im;

            if (dz_re * dz_re + dz_im * dz_im < 1e-12) {
                double r = sqrt(z_re * z_re + z_im * z_im);
                double theta = atan2(z_im, z_re);
                if (r > 1e-6) {
                    frequencies[found] = theta * fs / (2.0 * M_PI);
                    if (frequencies[found] < 0) frequencies[found] += fs;
                    damping[found] = log(r) * fs;
                    found++;
                }
                break;
            }
        }
    }
    return found;
}

#endif // PRONY_H
```

---

## 6. Подготовка данных для модального анализа

Hankel-матрица:

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

/**
 * @brief Построение Hankel-матрицы из временного ряда.
 *
 * @param signal    Входной сигнал длины N.
 * @param N         Длина сигнала.
 * @param num_rows  Число строк матрицы (1 <= num_rows <= N).
 * @param hankel    Выходной буфер размера num_rows * num_cols.
 * @return          Число столбцов матрицы (num_cols = N - num_rows + 1), или 0 при ошибке.
 *
 * @note Матрица хранится по строкам: H[i][j] = hankel[i * num_cols + j].
 */
size_t build_hankel_matrix(const double* signal, size_t N, size_t num_rows, double* hankel) {
    if (num_rows == 0 || num_rows > N) return 0;
    size_t num_cols = N - num_rows + 1;
    for (size_t i = 0; i < num_rows; ++i) {
        memcpy(&hankel[i * num_cols], &signal[i], num_cols * sizeof(double));
    }
    return num_cols;
}

/**
 * @brief Центрирование сигнала (удаление среднего значения).
 *
 * @param data Массив данных.
 * @param n    Длина массива.
 *
 * @note Модифицирует исходный массив in-place.
 */
void demean_signal(double* data, size_t n) {
    double sum = 0.0;
    for (size_t i = 0; i < n; ++i) sum += data[i];
    double mean = sum / (double)n;
    for (size_t i = 0; i < n; ++i) data[i] -= mean;
}

#endif // HANKEL_H
```

---

## 7. Использование

Пример полной обработки одного канала:

```c
#include "calibration.h"
#include "filtering.h"
#include "biquad.h"
#include "fht.h"
#include "prony.h"
#include "hankel.h"

#define FS 50.0
#define N 512
#define P 4

void full_processing_pipeline(
    const double* raw,
    const axis_calib_t* calib,
    double* spectrum_fht,
    double* freqs,
    double* damp
) {
    double x[N];
    for (int i = 0; i < N; ++i) {
        x[i] = apply_calibration(raw[i], calib);
    }

    hp_filter_state_t hp;
    hp_filter_init(&hp);
    double alpha = FS / (FS + 0.01);
    for (int i = 0; i < N; ++i) {
        x[i] = hp_filter_step(x[i], &hp, alpha);
    }

    biquad_filter_t bp;
    biquad_init(&bp, BIQUAD_BPF, 5.0, 5.0, FS);
    for (int i = 0; i < N; ++i) {
        x[i] = biquad_process(&bp, x[i]);
    }

    memcpy(spectrum_fht, x, N * sizeof(double));
    fht(spectrum_fht, N);

    int modes = prony_full(x, N, P, FS, freqs, damp);
    (void)modes;
}
```