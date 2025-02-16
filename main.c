#include <stdlib.h>
#include <stdio.h>

#include <math.h>
#include <complex.h>

double FFT_1D(double complex* w, int n, double k, double d, double theta) {
/* 
 * Given a 1D linear array of complex weights w with size n,
 * return the array response at the provided theta angle
 *
 * Args:
 * w - 1D array of complex element weights 
 * n - Number of elements in a 1D array slice
 * k - The wavenumber of the signal of interest
 * d - The distance in meters between elements
 * theta - The theta angle to evaluate at (-pi/2, pi/2)
 */

    double sum = 0;

    for (int i=0; i<n; i++) {
        sum += w[i] * cexp(-1*I * k * d * sin(theta) * i);
    }

    return sum;
}

double wavenumber(double f) {
/* Return the free space wavenumber given a frequency of f Hz */

    const double c = 299792458;

    double lambda = c/f;
    return (2 * M_PI) / lambda;
}

int main(void) {
    double complex* w = malloc(8 * sizeof(double complex));

    for (int i=0; i<8; i++) {
        w[i] = 1.0 + 0.0 * I; 

        printf("%.1f%+.1fi,", creal(w[i]), cimag(w[i]));
    }
    printf("\n");

    /* example: 55mm linear array at 2 GHz */
    double amp;

    for (double theta=0; theta < M_PI/2; theta += 0.1) {
        amp = cabs(FFT_1D(w, 8, wavenumber(2e9), 55e-3, theta));

        printf("theta=%.2f = %.2f\n", theta, amp);
    }

    free(w);

    return 0;
}
