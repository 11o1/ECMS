Ниже представлен **единый технический документ** уровня научно-инженерной разработки для системы мониторинга инженерных сооружений (мостов) на основе данных акселерометров и гироскопов. Документ включает:

- теоретические основы,
- математические формулы,
- высокоэффективные алгоритмы,
- реализацию на **чистом C99** (без внешних зависимостей, кроме `<math.h>`),
- комментарии, соответствующие стандартам профессиональной разработки.

---

# 📄 **Технический документ: Алгоритмы предварительной и спектральной обработки сигналов для систем мониторинга инженерных сооружений**

**Автор**: Система мониторинга AXIOMICA  
**Дата**: 23 октября 2025 г.  
**Цель**: Универсальная, вычислительно эффективная обработка данных от разнородных инерциальных датчиков с последующим модальным анализом.

---

## 1. Общая архитектура обработки

Сигнал проходит следующие этапы:

1. **Калибровка и компенсация смещения**  
2. **Фильтрация (подавление шума и дрейфа)**  
3. **Сегментация и оконное взвешивание**  
4. **Спектральное и параметрическое преобразование**  
5. **Подготовка данных для модального анализа (SSI/ESPRIT)**

Все этапы реализованы на C99 с учётом ограничений встраиваемых систем.

---

## 2. Калибровка датчиков

### 2.1. Модель ошибок MEMS-датчиков

Для акселерометра/гироскопа:

\[
x_{\text{raw}} = s \cdot x_{\text{true}} + b + n
\]

где:  
- \( s \) — масштабный коэффициент,  
- \( b \) — смещение нуля (bias),  
- \( n \) — аддитивный шум.

### 2.2. Реализация

```c
// calibration.c
#include "calibration.h"

double apply_calibration(double raw, const axis_calib_t* calib) {
    // y = (x - b) * s
    return (raw - calib->offset) * calib->scale;
}
```

> **Примечание**: Калибровочные коэффициенты определяются в лаборатории (метод наименьших квадратов по статическим измерениям).

---

## 3. Фильтрация сигналов

### 3.1. Высокочастотный фильтр (HPF) — подавление дрейфа

Передаточная функция аналогового HPF первого порядка:

\[
H(s) = \frac{s}{s + \omega_c}, \quad \omega_c = \frac{1}{\tau}
\]

Дискретизация по методу билинейного преобразования даёт:

\[
y[n] = \alpha \cdot (y[n-1] + x[n] - x[n-1]), \quad \alpha = \frac{2f_s - \omega_c}{2f_s + \omega_c}
\]

Однако для упрощения и устойчивости в условиях медленного дрейфа используется **интегрирующий HPF**:

\[
y[n] = \alpha \cdot (y[n-1] + x[n] - y[n-1])
\]

#### Реализация

```c
// filtering.c
#include "filtering.h"
#include <math.h>

void hp_filter_init(hp_filter_state_t* state) {
    state->y_prev = 0.0;
}

double hp_filter_step(double input, hp_filter_state_t* state, double alpha) {
    double y = alpha * (state->y_prev + input - state->y_prev);
    state->y_prev = y;
    return y;
}
```

> **Рекомендация**: Для мостов \( f_c = 0.01\,\text{Гц} \), при \( f_s = 50\,\text{Гц} \) → \( \alpha \approx 0.998 \).

---

### 3.2. Полосовой фильтр Баттерворта (Butterworth Bandpass)

Для выделения диапазона модальных частот (0.1–10 Гц).

Передаточная функция аналогового фильтра 2-го порядка:

\[
H(s) = \frac{s \cdot \omega_0 / Q}{s^2 + s \cdot \omega_0 / Q + \omega_0^2}
\]

Дискретизация — билинейное преобразование с предварительной коррекцией частоты:

\[
\omega_d = \tan\left( \frac{\pi f}{f_s} \right)
\]

#### Реализация IIR-фильтра 2-го порядка (биквад)

```c
typedef struct {
    double x1, x2;  // задержки входа
    double y1, y2;  // задержки выхода
    double b0, b1, b2; // числитель
    double a1, a2;     // знаменатель (a0 = 1)
} biquad_filter_t;

void biquad_init_lowpass(biquad_filter_t* f, double fc, double fs) {
    double omega = M_PI * fc / fs;
    double sn = sin(omega);
    double cs = cos(omega);
    double alpha = sn / (2.0 * 0.7071); // Q = 1/sqrt(2) для Баттерворта

    double b0 = (1.0 - cs) * 0.5;
    double b1 = 1.0 - cs;
    double b2 = b0;
    double a0 = 1.0 + alpha;
    double a1 = -2.0 * cs;
    double a2 = 1.0 - alpha;

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
```

> Для полосового фильтра комбинируются HPF и LPF или используется прямая формула полосового биквада.

---

## 4. Спектральные преобразования

### 4.1. Дискретное преобразование Хартли (DHT)

Определение:

\[
H[k] = \sum_{n=0}^{N-1} x[n] \cdot \operatorname{cas}\left( \frac{2\pi kn}{N} \right), \quad \operatorname{cas}(\theta) = \cos\theta + \sin\theta
\]

Преимущество: полностью вещественное, обратимое, энергосберегающее.

#### Реализация

```c
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
```

> **Ограничение**: \( O(N^2) \). Для \( N > 512 \) использовать FHT или FFT.

---

### 4.2. Быстрое преобразование Фурье (FFT) — Radix-2, in-place

Хотя DHT теоретически привлекателен, **FFT** предпочтителен из-за зрелости и оптимизации.

#### Реализация (базовая, для чётных N)

```c
// bit-reversal permutation
static void bit_reverse(double* x, size_t N) {
    for (size_t i = 0, j = 0; i < N; ++i) {
        if (j > i) {
            double tmp_real = x[2*i];   x[2*i] = x[2*j];   x[2*j] = tmp_real;
            double tmp_imag = x[2*i+1]; x[2*i+1] = x[2*j+1]; x[2*j+1] = tmp_imag;
        }
        size_t m = N >> 1;
        while (m >= 1 && j >= m) { j -= m; m >>= 1; }
        j += m;
    }
}

void fft_radix2(double* x, size_t N) {
    // x: массив длины 2*N: [Re0, Im0, Re1, Im1, ...]
    bit_reverse(x, N);
    for (size_t L = 2; L <= N; L <<= 1) {
        size_t L2 = L >> 1;
        double theta = -2.0 * M_PI / (double)L;
        double ur = 1.0, ui = 0.0;
        double wr = cos(theta), wi = sin(theta);
        for (size_t j = 0; j < L2; ++j) {
            for (size_t i = j; i < N; i += L) {
                size_t i1 = i + L2;
                double t_r = ur * x[2*i1] - ui * x[2*i1 + 1];
                double t_i = ur * x[2*i1 + 1] + ui * x[2*i1];
                x[2*i1]     = x[2*i] - t_r;
                x[2*i1 + 1] = x[2*i + 1] - t_i;
                x[2*i]     += t_r;
                x[2*i + 1] += t_i;
            }
            double tmp = ur * wr - ui * wi;
            ui = ur * wi + ui * wr;
            ur = tmp;
        }
    }
}
```

> **Примечание**: Для вещественных сигналов можно использовать оптимизацию (real FFT), но для общности приведена комплексная версия.

---

## 5. Параметрический анализ: метод Прони (Prony)

### 5.1. Математическая модель

Сигнал моделируется как сумма затухающих комплексных экспонент:

\[
x[n] = \sum_{k=1}^{p} A_k e^{(\sigma_k + j\omega_k) n T}
\]

Цель — найти \( \omega_k \) (собственные частоты) и \( \sigma_k \) (декременты).

### 5.2. Алгоритм

1. Решить систему \( \mathbf{X} \mathbf{a} = -\mathbf{x} \), где \( \mathbf{X} \) — автокорреляционная матрица.
2. Найти корни полинома \( P(z) = 1 + a_1 z^{-1} + \dots + a_p z^{-p} \).
3. Извлечь частоты: \( \omega_k = \arg(z_k) / T \).

#### Реализация (упрощённая, для малого p)

```c
// prony.c
#include <math.h>
#include <string.h>

// Решение системы методом наименьших квадратов (нормальные уравнения)
static void solve_least_squares(const double* X, const double* y, double* a, size_t rows, size_t cols) {
    // X: rows x cols, y: rows, a: cols
    // Формируем A = X^T X, b = X^T y
    double A[cols * cols];
    double b[cols];
    memset(A, 0, sizeof(A));
    memset(b, 0, sizeof(b));

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            b[j] += X[i * cols + j] * y[i];
            for (size_t k = 0; k < cols; ++k) {
                A[j * cols + k] += X[i * cols + j] * X[i * cols + k];
            }
        }
    }

    // Решаем A a = b методом Гаусса (упрощённо)
    // ... (реализация Гаусса опущена для краткости; в production — использовать QR или SVD)
}

int prony_estimate_frequencies(const double* x, size_t N, size_t p,
                               double fs, double* frequencies, double* damping) {
    if (p >= N) return -1;

    // Формируем матрицу X (размер (N-p) x p)
    double* X = (double*)malloc((N - p) * p * sizeof(double));
    double* y = (double*)malloc((N - p) * sizeof(double));
    for (size_t i = 0; i < N - p; ++i) {
        for (size_t j = 0; j < p; ++j) {
            X[i * p + j] = x[i + p - 1 - j]; // обратный порядок
        }
        y[i] = -x[i + p];
    }

    double* a = (double*)calloc(p, sizeof(double));
    solve_least_squares(X, y, a, N - p, p);

    // Теперь находим корни полинома: z^p + a[0] z^{p-1} + ... + a[p-1] = 0
    // Это требует реализации поиска корней полинома (например, метод Дженкинса-Трауба)
    // В рамках документа — заглушка:
    // Предположим, что корни уже найдены как complex_roots[0..p-1]

    // Пример для одного корня (демонстрация):
    // double r = sqrt(creal*creal + cimag*cimag);
    // double theta = atan2(cimag, creal);
    // frequencies[0] = theta * fs / (2*M_PI);
    // damping[0] = log(r) * fs;

    free(X); free(y); free(a);
    return 0; // успех (в реальной реализации — число найденных мод)
}
```

> **Рекомендация**: Для надёжности использовать **TLS-ESPRIT** или **Stochastic Subspace Identification (SSI)** на сервере.

---

## 6. Подготовка данных для SSI (Stochastic Subspace Identification)

### 6.1. Построение Hankel-матрицы

Для сигнала \( y[0], \dots, y[N-1] \):

\[
\mathbf{H} = \begin{bmatrix}
y[0] & y[1] & \cdots & y[j-1] \\
y[1] & y[2] & \cdots & y[j] \\
\vdots & \vdots & \ddots & \vdots \\
y[i-1] & y[i] & \cdots & y[i+j-2]
\end{bmatrix}
\]

где \( i + j - 1 = N \).

#### Реализация

```c
size_t build_hankel_matrix(const double* signal, size_t N, size_t num_rows, double* hankel) {
    if (num_rows == 0 || num_rows > N) return 0;
    size_t num_cols = N - num_rows + 1;
    for (size_t i = 0; i < num_rows; ++i) {
        memcpy(&hankel[i * num_cols], &signal[i], num_cols * sizeof(double));
    }
    return num_cols;
}
```

> Эта матрица подаётся на вход SVD в рамках SSI.

---

## 7. Заключение и рекомендации

| Этап | Метод | Применение |
|------|-------|-----------|
| Калибровка | Аффинная модель | Обязательно |
| Фильтрация | HPF + Butterworth BP | Подавление дрейфа и шума |
| Спектр | FFT / DHT | Быстрая оценка спектра |
| Модальный анализ | Prony / SSI / ESPRIT | Точный выделение собственных частот |
| Трекинг | Фильтр Калмана | Мониторинг изменений во времени |

### Требования к реализации:
- На **контроллере**: HPF, калибровка, буферизация.
- На **сервере**: SSI, кластеризация, визуализация.

### Дальнейшие шаги:
1. Реализовать **SVD на чистом C** (Golub-Reinsch).
2. Интегрировать **мультидатчиковую синхронизацию**.
3. Добавить **температурную компенсацию** модальных частот.

---

**Приложение**: Все заголовочные файлы (`*.h`) и реализации (`*.c`) совместимы с **C99**, не используют динамическое выделение памяти (кроме `prony.c` — его можно переписать под статические буферы), и могут быть скомпилированы с `-O2 -ffast-math` для встраиваемых систем.

Если требуется полная реализация **SVD**, **ESPRIT** или **серверного конвейера**, сообщите — предоставлю в том же формате.
