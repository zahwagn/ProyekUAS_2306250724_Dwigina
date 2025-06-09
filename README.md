# Aplikasi Program Range-Kutta Orde ke-4 pada Dekomposisi Kimia
```text
Nama    : Dwigina Sitti Zahwa
NPM     : 2306250724
```

## Deskripsi Program
Program ini mensimulasikan dekomposisi senyawa kimia yang memiliki **konsentrasi awal 10 juta partikel/m³** menggunakan metode Runge-Kutta Orde 4 (RK4) untuk menghitung **perubahan konsentrasi setelah 7 minggu** dengan **konstanta laju reaksi 0.2 per minggu**. Program ini menghitung dengan RK4 dengan persamaan diferensial orde pertama:
$$\frac{dC}{dt} = -kC$$

- $C$: konsentrasi awal (parts/m³)
- $t$: waktu (minggu)
- $k$: konstanta laju reaksi (0.2 per minggu)

Program memberikan *output* antara lain:
- Prediksi konsentrasi setiap minggu sampai minggu ke-7 menggunakan RK4
- Perbandingan dengan solusi analitik -> $C(t) = C_0 e^{-kt}$
- Perhitungan error numerik (True Error dan Approximate Error)

## Metode 
Metode numerik **Runge-Kutta Orde 4 (RK4)** merupakan metode yang biasanya dipakai untuk menyelesaikan masalah nilai awal (initial value problem) dari persamaan diferensial biasa orde satu. Metode ini memberikan akurasi tinggi dengan kesalahan global sebesar $O(h^5)$. Metode RK4 akan melalui evaluasi empat ***slope*** ($k1$ - $k4$) yang didefinisikan:

$$k_1 = f(t_n, y_n) \\,k_2 = f(t_n + \tfrac{h}{2}, y_n + \tfrac{hk_1}{2}) \\,k_3 = f(t_n + \tfrac{h}{2}, y_n + \tfrac{hk_2}{2}) \\,k_4 = f(t_n + h, y_n + hk_3)$$

Lalu update:
$$y_{n+1} = y_n + \frac{h}{6}(k_1 + 2k_2 + 2k_3 + k_4)$$

## Pseudocode
#### (a) Main or "Driver" Program

```text
Assign values for:
n = number of equations
yi = initial values of n dependent variables
xi = initial value independent variable
xf = final value independent variable
dx = calculation step size
xout = output interval

x = xi
m = 0
xpm = x

DOFOR i = 1 TO n
    ypi,m = yii
    yi = yii
END DO

DO
    xend = x + xout
    IF (xend > xf) THEN xend = xf
    h = dx
    CALL Integrator(x, y, n, h, xend)
    m = m + 1
    xpm = x

    DOFOR i = 1 TO n
        ypi,m = yi
    END DO

    IF (x >= xf) EXIT
END DO

DISPLAY RESULTS
```

#### (b) Routine to Take One Output Step

```text
SUBROUTINE Integrator(x, y, n, h, xend)
    DO
        IF (xend - x < h) THEN
            h = xend - x
        END IF
        CALL RK4(x, y, n, h)
        IF (x >= xend) EXIT
    END DO
END SUBROUTINE
```
#### (c) Fourth-Order RK Method for a System of ODEs

```text
SUBROUTINE RK4(x, y, n, h)
    CALL Derivs(x, y, k1)

    DOFOR i = 1 TO n
        ymi = yi + k1i * h / 2
    END DO

    CALL Derivs(x + h/2, ym, k2)

    DOFOR i = 1 TO n
        ymi = yi + k2i * h / 2
    END DO

    CALL Derivs(x + h/2, ym, k3)

    DOFOR i = 1 TO n
        yei = yi + k3i * h
    END DO

    CALL Derivs(x + h, ye, k4)

    DOFOR i = 1 TO n
        slopei = (k1i + 2*(k2i + k3i) + k4i)/6
        yi = yi + slopei * h
    END DO

    x = x + h
END SUBROUTINE
```

#### (d) Routine to Determine Derivatives

```text
SUBROUTINE Derivs(x, y, dy)
    dy1 = ...   ! Fungsi turunan untuk variabel pertama
    dy2 = ...   ! Fungsi turunan untuk variabel kedua (jika ada)
END SUBROUTINE
```

## Source Code
```c 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Jumlah persamaan diferensial 
#define NE 1

/**
 * Fungsi turunan untuk persamaan diferensial
 * x = Variabel independen (waktu)
 * y = Array variabel dependen (konsentrasi saat ini)
 * dydx =  Array output turunan
 */
void Derivs(double x, double y[], double dydx[]) {
    double k = 0.2; // Konstanta laju reaksi (per minggu)
    dydx[0] = -k * y[0]; // Persamaan diferensial orde pertama = dC/dt
}

/**
 * Implementasi metode Runge-Kutta Orde 4
 * x = waktu saat ini
 * y = Nkonsentrasi saat ini
 * h = Step-size
 * n = Jumlah persamaan (1)
 * Derivs = Fungsi turunan
 */
void RK4(double x, double y[], double h, int n, void (*Derivs)(double, double[], double[])) {
    double k1[NE], k2[NE], k3[NE], k4[NE]; // Kemiringan (slope) untuk RK4
    double ytemp[NE]; // Nilai sementara variabel dependen

    // Tahapan-1 evaluasi slope  k1
    Derivs(x, y, k1);
    for (int i = 0; i < NE; i++) {
        ytemp[i] = y[i] + 0.5 * h * k1[i];
    }

// Tahapan-2 evaluasi slope  k2
    Derivs(x + 0.5 * h, ytemp, k2);
    for (int i = 0; i < NE; i++) {
        ytemp[i] = y[i] + 0.5 * h * k2[i];
    }

    // Tahapan-3 evaluasi slope  k3
    Derivs(x + 0.5 * h, ytemp, k3);
    for (int i = 0; i < NE; i++) {
        ytemp[i] = y[i] + h * k3[i];
    }

    // Tahapan-4 evaluasi slope  k4
    Derivs(x + h, ytemp, k4);

    // Gabungkan semua slope untuk update nilai y
    for (int i = 0; i < NE; i++) {
        y[i] = y[i] + h*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i])/6.0;
    }
}

/**
 * Integrator untuk mengelola step-size waktu
 * x = Pointer ke waktu saat ini (update setelah setiap langkah)
 * y =  Konsentrasi yang akan diupdate
 * h = step-size awal
 * xend = Waktu akhir untuk integrasi saat ini
 * n = Jumlah persamaan
 * RK4 = Fungsi metode RK4
 */
void Integrator(double *x, double y[], double h, double xend, int n, 
                void (*RK4)(double, double[], double, int, void (*)(double, double[], double[]))) {
    while ((*x) < xend) {
        if (xend - *x < h) h = xend - *x; // Sesuaikan step-size di akhir
        RK4(*x, y, h, n, Derivs); // Panggil RK4
        (*x) += h; // Update waktu
    }
}

int main() {
    // Inisialisasi parameter
    double xi = 0.0;          // Waktu awal (minggu)
    double xf = 7.0;          // Waktu akhir (minggu)
    double dx = 1.0;          // step-size (minggu)
    double xout = 1.0;        // Interval output
    double C0 = 1e7;          // Konsentrasi awal (partikel/m³)
    double y[NE] = {C0};      // Variabel dependen
    double x = xi;            // Variabel independen
    int m = 0;                // Counter output
    double xpm[100];          // Penyimpanan waktu output
    double ypm[NE][100];      // Penyimpanan konsentrasi output

    // Simpan kondisi awal
    xpm[m] = x;
    for (int i = 0; i < NE; i++) ypm[i][m] = y[i];

    // Loop integrasi utama
    while (x < xf) {
        double xend = x + xout;
        if (xend > xf) xend = xf;
        Integrator(&x, y, dx, xend, NE, RK4);
        m++;
        xpm[m] = x;
        for (int i = 0; i < NE; i++) ypm[i][m] = y[i];
    }

    // Output hasil dan error
    printf("Minggu\tKonsentrasi (RK4)\tNilai Eksak\tTrue Error (%)\tApprox Error (%)\n");
    for (int i = 0; i <= m; i++) {
        double true_val = C0 * exp(-0.2 * xpm[i]); // Solusi analitik
        double true_error = fabs((true_val - ypm[0][i]) / true_val) * 100;
        double approx_error = (i == 0) ? 0.0 : fabs((ypm[0][i] - ypm[0][i-1]) / ypm[0][i]) * 100;
        
        printf("%.1f\t%.2f\t\t%.2f\t\t%.4f\t\t%.4f\n", 
               xpm[i], ypm[0][i], true_val, true_error, approx_error);
    }

    return 0;
}
```
## Penjelasan Proses

#### **1. Fungsi `Derivs` (Subrutin Turunan)**
**Pseudocode**:  
```plaintext
SUB Derivs(x, y, dy)
    dy1 = -k*y1
END SUB
```

**Implementasi & Penjelasan**:
```c
void Derivs(double x, double y[], double dydx[]) {
    double k = 0.2; // Konstanta laju reaksi
    dydx[0] = -k * y[0]; // Menghitung turunan konsentrasi (dC/dt = -kC)
}
```
- **Tujuan**: Menghitung turunan konsentrasi berdasarkan persamaan diferensial orde pertama.
- **Proses**:
  - Menghitung `dC/dt` menggunakan hukum laju reaksi kimia orde pertama (`-k*C`).

---

#### **2. Fungsi `RK4` (Metode Runge-Kutta Orde 4)**
**Pseudocode**:  
```plaintext
SUB RK4(x, y, n, h)
    CALL Derivs(x, y, k1)
    DOFOR i = 1, n
        ymi = yi + k1i * h / 2
    END DO
    CALL Derivs(x + h/2, ym, k2)
    ... (Langkah k3 dan k4 serupa)
    DOFOR i = 1, n
        slopei = (k1i + 2*(k2i + k3i) + k4i)/6
        yi = yi + slopei * h
    END DO
END SUB
```

**Implementasi & Penjelasan**:
```c
void RK4(double x, double y[], double h, int n, void (*Derivs)(double, double[], double[])) {
    double k1[NE], k2[NE], k3[NE], k4[NE]; // Kemiringan (slope)
    double ytemp[NE]; // Nilai sementara
```
- **Tujuan**: Mengimplementasikan metode Runge-Kutta Orde 4 dengan evaluasi 4 titik slope.
- **Proses**:
  1. **k1**: Kemiringan di titik awal.
  2. **k2**: Kemiringan di titik tengah menggunakan `k1`.
  3. **k3**: Kemiringan di titik tengah menggunakan `k2`.
  4. **k4**: Kemiringan di titik akhir menggunakan `k3`.
  5. **Update `y`**: nilai rata-rata dari keempat slope.

---

#### **3. Fungsi `Integrator` (Pengelola Langkah Waktu)**
**Pseudocode**:  
```plaintext
SUB Integrator(x, y, n, h, xend)
    DO
        IF (xend - x < h) THEN h = xend - x
        CALL RK4(x, y, n, h)
        IF (x >= xend) EXIT
    END DO
END SUB
```

**Implementasi & Penjelasan**:
```c
void Integrator(double *x, double y[], double h, double xend, int n, 
                void (*RK4)(double, double[], double, int, void (*)(double, double[], double[]))) {
    while (*x < xend) {
        if (xend - *x < h) h = xend - *x; // Sesuaikan langkah terakhir
        RK4(*x, y, h, n, Derivs); // Panggil RK4
        *x += h; // Update waktu
    }
}
```
- **Tujuan**: Mengintegrasikan persamaan dari waktu `x` ke `xend` dengan step `h`.
- **Proses**:
  - Loop hingga mencapai `xend`.
  - Menyesuaikan step-size jika mendekati `xend`.
  - Memanggil `RK4` untuk setiap langkah.

---

#### **4. Fungsi `main` (Driver Program)**
**Pseudocode**:  
```plaintext
MAIN PROGRAM
    Assign nilai awal: xi, xf, dx, C0
    DO WHILE x < xf
        xend = x + xout
        CALL Integrator(x, y, n, dx, xend)
        Simpan hasil
    END DO
    Tampilkan hasil
END PROGRAM
```

**Implementasi & Penjelasan**:
```c
int main() {
    // Inisialisasi parameter
    double xi = 0.0, xf = 7.0, dx = 1.0, xout = 1.0;
    double C0 = 1e7, y[NE] = {C0}, x = xi;
```
- **Tujuan**: Mengkoordinasi seluruh proses simulasi.
- **Proses**:
  1. Inisialisasi waktu, konsentrasi awal, dan parameter lain.
  2. Loop integrasi dengan memanggil `Integrator` untuk setiap interval output.
  3. Menyimpan dan menampilkan hasil perbandingan dengan solusi eksak.

---
