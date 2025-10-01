#include "signal_processing.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void write_array_to_file(const char *filename, const double *data, int size) {
    FILE *f = fopen(filename, "w");
    if (!f) {
        perror("Failed to open file for writing");
        return;
    }
    for (int i = 0; i < size; i++) {
        fprintf(f, "%.8f\n", data[i]);
    }
    fclose(f);
}

void calculate_autocorrelation(const SignalData *signal, double *acf_result) {
    int n = signal->count;
    for (int lag = 0; lag < n; lag++) {
        double sum = 0.0;
        for (int i = 0; i < n - lag; i++) {
            sum += signal->data[i] * signal->data[i + lag];
        }
        acf_result[lag] = sum / n;
    }
}

void psd_periodogram(const SignalData *signal, double *psd, fftw_plan *p) {
    int n = signal->count;
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (n / 2 + 1));
    if (*p == NULL) {
        *p = fftw_plan_dft_r2c_1d(n, (double*)signal->data, out, FFTW_ESTIMATE);
    }
    fftw_execute(*p);
    for (int i = 0; i < n / 2 + 1; i++) {
        psd[i] = (out[i][0] * out[i][0] + out[i][1] * out[i][1]) / n;
    }
    fftw_free(out);
}

void psd_welch(const SignalData *signal, double *psd, int nperseg, int noverlap) {
    int n = signal->count;
    int step = nperseg - noverlap;
    if (step <= 0) return;
    int num_segments = (n - noverlap) / step;
    if (num_segments == 0) { return; }

    double *segment = (double*) fftw_malloc(sizeof(double) * nperseg);
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (nperseg / 2 + 1));
    double *window = (double*) malloc(sizeof(double) * nperseg);
    for(int i = 0; i < nperseg; i++) {
        window[i] = 0.5 * (1 - cos(2 * M_PI * i / (nperseg - 1))); // Hanning window
    }

    fftw_plan p = fftw_plan_dft_r2c_1d(nperseg, segment, out, FFTW_ESTIMATE);
    memset(psd, 0, sizeof(double) * (nperseg / 2 + 1));

    for (int i = 0; i < num_segments; i++) {
        for (int j = 0; j < nperseg; j++) {
            segment[j] = signal->data[i * step + j] * window[j];
        }
        fftw_execute(p);
        for (int j = 0; j < nperseg / 2 + 1; j++) {
            psd[j] += (out[j][0] * out[j][0] + out[j][1] * out[j][1]) / nperseg;
        }
    }
    for (int i = 0; i < nperseg / 2 + 1; i++) {
        psd[i] /= num_segments;
    }

    fftw_destroy_plan(p);
    fftw_free(segment);
    fftw_free(out);
    free(window);
}

void psd_multitaper(const SignalData *signal, double *psd) {
    int n = signal->count;
    int num_tapers = 4;
    double *tapered_signal = (double*) fftw_malloc(sizeof(double) * n);
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (n / 2 + 1));
    fftw_plan p = fftw_plan_dft_r2c_1d(n, tapered_signal, out, FFTW_ESTIMATE);
    memset(psd, 0, sizeof(double) * (n / 2 + 1));

    for (int k = 0; k < num_tapers; k++) {
        for (int i = 0; i < n; i++) {
            tapered_signal[i] = signal->data[i] * sin(M_PI * (i + 0.5) / n);
        }
        fftw_execute(p);
        for (int i = 0; i < n / 2 + 1; i++) {
            psd[i] += (out[i][0] * out[i][0] + out[i][1] * out[i][1]) / n;
        }
    }
    for (int i = 0; i < n / 2 + 1; i++) {
        psd[i] /= num_tapers;
    }

    fftw_destroy_plan(p);
    fftw_free(tapered_signal);
    fftw_free(out);
}