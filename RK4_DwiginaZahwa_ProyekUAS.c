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
