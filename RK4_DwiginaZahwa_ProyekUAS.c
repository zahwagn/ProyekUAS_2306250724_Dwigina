#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NE 1

// Differential equation: dC/dt = -kC
void Derivs(double x, double y[], double dydx[]) {
    double k = 0.2; // Degradation rate constant (per week)
    dydx[0] = -k * y[0]; // First-order degradation
}

// 4th Runge-Kutta Method
void RK4(double x, double y[], double h, int n, void (*Derivs)(double, double[], double[])) {
    double k1[NE], k2[NE], k3[NE], k4[NE];
    double ytemp[NE];

    // Calculate k1
    Derivs(x, y, k1);
    int i;
    for (i = 0; i < NE; i++) {
        ytemp[i] = y[i] + 0.5 * h * k1[i];
    }

    // Calculate k2
    Derivs(x + 0.5 * h, ytemp, k2);
    for (i = 0; i < NE; i++) {
        ytemp[i] = y[i] + 0.5 * h * k2[i];
    }

    // Calculate k3
    Derivs(x + 0.5 * h, ytemp, k3);
    for (i = 0; i < NE; i++) {
        ytemp[i] = y[i] + h * k3[i];
    }

    // Calculate k4
    Derivs(x + h, ytemp, k4);

    // Update y and x
    for (i = 0; i < NE; i++) {
        y[i] = y[i] + h*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i])/6.0;
    }
}

// Integrator to take one output step
void Integrator(double *x, double y[], double h, double xend, int n, 
                void (*RK4)(double, double[], double, int, void (*)(double, double[], double[]))) {
    while ((*x) < xend) {
        if (xend - *x < h) h = xend - *x; // Adjust step size 
        RK4(*x, y, h, n, Derivs);
        (*x) += h;
    }
}

int main() {
    // Parameters
    double xi = 0.0;          // Initial time (weeks)
    double xf = 7.0;          // Final time (weeks)
    double dx = 1.0;          // Step size (weeks)
    double xout = 1.0;        // Output interval
    double C0 = 1e7;          // Initial concentration
    double y[NE] = {C0};      // Dependent variables
    double x = xi;
    int m = 0;                
    double xpm[100];          
    double ypm[NE][100];      

    // Save initial condition
    xpm[m] = x;
    int i;
    for (i = 0; i < NE; i++) ypm[i][m] = y[i];

    // Main integration loop
    while (x < xf) {
        double xend = x + xout;
        if (xend > xf) xend = xf;
        Integrator(&x, y, dx, xend, NE, RK4);
        m++;
        xpm[m] = x;
        for (i = 0; i < NE; i++) ypm[i][m] = y[i];
    }

    // Output results and errors
    printf("Week\tConcentration (RK4)\tTrue Value\tTrue Error (%)\tApproximate Error (%)\n");
    for (i = 0; i <= m; i++) {
        double true_val = C0 * exp(-0.2 * xpm[i]);
        double approx_error = (i == 0) ? 0.0 : fabs((ypm[0][i] - ypm[0][i-1]) / ypm[0][i]) * 100;
        double true_error = fabs((true_val - ypm[0][i]) / true_val) * 100;
        printf("%.1f\t%.2f\t\t%.2f\t\t%.4f\t\t%.4f\n", 
               xpm[i], ypm[0][i], true_val, true_error, approx_error);
    }

    return 0;
}
