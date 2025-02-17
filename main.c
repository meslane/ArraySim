#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <math.h>
#include <complex.h>

struct phased_array {
    double complex** w; //weights
    int m; //number of elements in x axis
    int n; //number of elements in y axis
    double k; //wavenumber
    double d; //element spacing
};


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


double FFT_2D(double complex** w, int n_x, int n_y, 
                double k, double d, 
                double el, double az) {
/* 
 * Given a 2D rectangular array of complex weights with size n,
 * return the array response at the provided theta and phi angle
 *
 * Args:
 * w - 2D array of complex element weights
 * n_x - Number of elements in the x array dimension
 * n_y - Number of elements in the y array dimension 
 * k - The wavenumber of the signal of interest
 * d - the distance in meters between elements (right now assumes uniform in x AND y)
 * el - The elevation angle to evaluate at
 * az - The azimuth angle to evaluate at
 */

    double sum = 0;

    double n_0 = ((double) n_x - 1.0) / 2.0;
    double m_0 = ((double) n_y - 1.0) / 2.0;

    for (int m=0; m < n_x; m++) {
        for (int n=0; n < n_y; n++) {
            sum += w[m][n] * cexp(I*k * d * ((m-n_0)*cos(az)*sin(el) + (n-m_0)*sin(az)*sin(el)));
        }
    }

    return sum;
}

double complex** alloc_array_2D(int m, int n, double complex val) {
/* Allocate a 2D complex array of size m * n 
 *
 * m - The size of the first array dimension
 * n - The size of the second array dimension
 * val - The value to initialize every element to
 */

    double complex** w = malloc(m * sizeof(double complex*));

    for (int i=0; i<m; i++) {
        w[i] = malloc(n * sizeof(double complex));
        for (int j=0; j<n; j++) {
            w[i][j] = val;
        }
    }

    return w;
}

void free_array_2D(double complex** w, int m) {
/* Free a 2D complex array w of size m * n
 *
 * w - The 2D array to be freed
 * m - The size of the first array dimension
 * n - The size of the second array dimension
 */

    for (int i=0; i<m; i++) {
        free(w[i]);
    }

    free(w);
}

void plot_pattern_cut(struct phased_array* a, double phi) {
/* Plot a 2D pattern cut at the given azimith angle
 *
 * a = Phased array to plot
 * phi = Phi angle at which to cut (in radians)
 */
    
    double scale = 1; 
    double norm = a->m * a->n;

    double amp;

    printf("%d\n", a->m);
    printf("%.1f%+.1fi\n", creal(a->w[1][1]), cimag(a->w[1][1]));

    for (int i = 0; i <= 30; i++) {
        for (double j = -1.5; j < 1.52; j += 0.02) {
            amp = 20*log10(cabs(FFT_2D(a->w, a->m, a->n, a->k, a->d, j, phi))/norm);          
            
            if (amp >= -1 * i / scale) {
                printf("*");
            }
            else {
                printf(" ");
            }
        }
        printf("\n");
    }
}

void plot_uv_pattern(struct phased_array* a, double range, double step) {
/* Plot the 2D phased array pattern in UV space
 *
 * a = Phased array to plot
 * range = Theta range to plot (max = 1.0)
 * step = UV step for each character drawn to the terminal
 */

    const static char amp_chars[][32] = {"\e[0;31m@\e[0m", "\e[0;33m%\e[0m", "\e[0;33m#\e[0m", "\e[0;93m*\e[0m", "\e[0;32m+\e[0m", "\e[0;32m=\e[0m", "\e[0;36m-\e[0m", "\e[0;34m:\e[0m", "\e[0;35m.\e[0m"}; 
    const static double amp_levels[] = {-3, -5, -10, -15, -20, -25, -30, -35, -40};

    double norm = a->m * a->n;
    double amp;
    double theta;
    double phi;
    static char draw_char[32] = {0};

    for (double u = -1*range; u < range+step; u += step) {
        for (double v = -1*range; v < range+step; v += step) {
            theta = asin(sqrt(u*u + v*v));
            phi = atan2(v,u);
            
            amp = 20*log10(cabs(FFT_2D(a->w, a->m, a->n, a->k, a->d, theta, phi))/norm);          
           
            //draw_char = " ";
            draw_char[0] = '\0';
            strncat(draw_char, " ", 30);
            for (int i = sizeof(amp_chars) / sizeof(char); i >= 0; i--) {
                if (amp >= amp_levels[i]) {
                    draw_char[0] = '\0';
                    strncat(draw_char, amp_chars[i], 30);
                }
            }

            printf("%s ", draw_char);
        }
        printf("\n");
    } 

}

struct phased_array* create_array(int m, int n, double freq, double d) {
/* Create a phased array struct
 *
 * m = Number of x elements
 * n = Number of y elements
 * freq = Frequency of interest (Hz)
 * d = Distance between elements (meters)
 */

    struct phased_array* a = malloc(sizeof(struct phased_array));

    a->m = m;
    a->n = n;
    a->w = alloc_array_2D(m, n, 1);
    a->k = wavenumber(freq);
    a->d = d;

    return a;
}

void free_array(struct phased_array* a) {
/* free memory allocated to phased array a */
    free_array_2D(a->w, a->m);
    free(a);
}

double complex ang_to_imag(double mag, double phase) {
/* Convert a magnitue and phase angle to a complex number 
 *
 * mag = magnitude
 * phase = phase (radians)
 */
    double r = mag * cos(phase);
    double im = mag * sin(phase);

    return CMPLX(r, im);
}

void steer_array(struct phased_array* a, double theta, double phi) {
/* Apply weights to steer the phased array to the given theta and phi angle
 *
 * a - Pointer to array to steer
 * theta - Theta pointing angle
 * phi - Phi pointing angle
 */

    double n_0 = ((double) a->n - 1.0) / 2.0;
    double m_0 = ((double) a->m - 1.0) / 2.0;

    double s = a->d * a->k;//steering scaling factor (d/lambda)

    for (int m=0; m < a->m; m++) {
        for (int n=0; n < a->n; n++) {
            a->w[m][n] = ang_to_imag(1.0, s * sin(theta) * ((cos(phi) * (m - m_0)) + (sin(phi) * (n - n_0))));
        }
    }
}

int main(int argc, char** argv) {
/* args:
 * 1: frequency (MHz)
 * 2: spacing
 * 3: x elements
 * 4: y elements
 * 5: plot resolution
 */

    if (argc < 7) {
        printf("\nCommand Line Arguments:\n");
        printf("1: frequency in Hz\n");
        printf("2: element spacing in meters\n");
        printf("3: x dimension element count\n");
        printf("4: y dimension element connt\n");
        printf("5: theta scan angle\n");
        printf("6: phi scan angle\n");
        printf("7: u/v resolution of the plot (0.05 for medium, 0.02 for bigger, etc.)\n\n");
        return 1;
    }

    struct phased_array* a = create_array(atof(argv[3]), atof(argv[4]), atof(argv[1]), atof(argv[2]));

    /*
    for (int i=0; i<3; i++) {
        for (int j = 0; j<3; j++) {
            a->w[i][j] = 0;
        }
    }
    */

    steer_array(a, atof(argv[5]), atof(argv[6]));

    plot_uv_pattern(a, 1.0, atof(argv[7]));

    free_array(a);

    return 0;
}
