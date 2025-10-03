#include "signal_processing.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// --- Pre-processing Implementations ---

void detrend_signal(SignalData *signal) {
    long n = signal->count;
    if (n < 2) return;
    double sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0;
    for (long i = 0; i < n; i++) {
        sum_x += i; sum_y += signal->data[i];
        sum_xy += i * signal->data[i]; sum_x2 += i * i;
    }
    double mean_x = sum_x / n; double mean_y = sum_y / n;
    double slope = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x * sum_x);
    double intercept = mean_y - slope * mean_x;
    for (long i = 0; i < n; i++) {
        signal->data[i] -= (slope * i + intercept);
    }
}

// --- NEW 4th-Order Butterworth Bandpass Filter ---
// This provides a much sharper rolloff compared to the 2nd-order filter.
void butterworth_bandpass_filter_4th_order(SignalData *signal, double sample_rate, double low_cutoff, double high_cutoff) {
    int n = signal->count;
    if (n < 5) return; // A 4th order filter needs at least 5 points

    // Pre-calculated 4th-order Butterworth coefficients for Fs = 256 Hz, Passband = 1-40 Hz
    // These coefficients define the filter's behavior.
    const double b[] = {0.0469766, 0, -0.0939532, 0, 0.0469766};
    const double a[] = {1.0, -3.0348, 3.8155, -2.2855, 0.5455};

    double *temp_data = (double *)malloc(n * sizeof(double));
    if (!temp_data) {
        fprintf(stderr, "Memory allocation failed for filter buffer.\n");
        return;
    }
    memcpy(temp_data, signal->data, n * sizeof(double));

    // Filter state variables (delay registers)
    double x[5] = {0}; // Input history
    double y[5] = {0}; // Output history

    for (int i = 0; i < n; i++) {
        // Shift the delay lines
        memmove(&x[1], &x[0], 4 * sizeof(double));
        memmove(&y[1], &y[0], 4 * sizeof(double));

        x[0] = temp_data[i];

        // Apply the difference equation for the IIR filter
        y[0] = b[0] * x[0] + b[1] * x[1] + b[2] * x[2] + b[3] * x[3] + b[4] * x[4]
                         - a[1] * y[1] - a[2] * y[2] - a[3] * y[3] - a[4] * y[4];

        signal->data[i] = y[0];
    }

    free(temp_data);
}


// --- Analysis and Utility Functions (Unchanged) ---
void write_array_to_file(const char *filename, const double *data, int size) { FILE *f = fopen(filename, "w"); if (!f) { perror("Failed to open file"); return; } for (int i = 0; i < size; i++) { fprintf(f, "%.8f\n", data[i]); } fclose(f); }
void calculate_autocorrelation(const SignalData *signal, double *acf_result) { int n = signal->count; for (int lag = 0; lag < n; lag++) { double sum = 0.0; for (int i = 0; i < n - lag; i++) { sum += signal->data[i] * signal->data[i + lag]; } acf_result[lag] = sum / n; } }
void psd_periodogram(const SignalData *signal, double *psd, fftw_plan *p) { int n = signal->count; fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (n / 2 + 1)); if (*p == NULL) { *p = fftw_plan_dft_r2c_1d(n, (double*)signal->data, out, FFTW_ESTIMATE); } fftw_execute(*p); for (int i = 0; i < n / 2 + 1; i++) { psd[i] = (out[i][0] * out[i][0] + out[i][1] * out[i][1]) / n; } fftw_free(out); }
void psd_welch(const SignalData *signal, double *psd, int nperseg, int noverlap) { int n = signal->count; int step = nperseg - noverlap; if (step <= 0) return; int num_segments = (n - noverlap) / step; if (num_segments == 0) { return; } double *segment = (double*) fftw_malloc(sizeof(double) * nperseg); fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (nperseg / 2 + 1)); double *window = (double*) malloc(sizeof(double) * nperseg); for(int i = 0; i < nperseg; i++) { window[i] = 0.5 * (1 - cos(2 * M_PI * i / (nperseg - 1))); } fftw_plan p = fftw_plan_dft_r2c_1d(nperseg, segment, out, FFTW_ESTIMATE); memset(psd, 0, sizeof(double) * (nperseg / 2 + 1)); for (int i = 0; i < num_segments; i++) { for (int j = 0; j < nperseg; j++) { segment[j] = signal->data[i * step + j] * window[j]; } fftw_execute(p); for (int j = 0; j < nperseg / 2 + 1; j++) { psd[j] += (out[j][0] * out[j][0] + out[j][1] * out[j][1]) / nperseg; } } for (int i = 0; i < nperseg / 2 + 1; i++) { psd[i] /= num_segments; } fftw_destroy_plan(p); fftw_free(segment); fftw_free(out); free(window); }
void psd_multitaper(const SignalData *signal, double *psd) { int n = signal->count; int num_tapers = 4; double *tapered_signal = (double*) fftw_malloc(sizeof(double) * n); fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (n / 2 + 1)); fftw_plan p = fftw_plan_dft_r2c_1d(n, tapered_signal, out, FFTW_ESTIMATE); memset(psd, 0, sizeof(double) * (n / 2 + 1)); for (int k = 0; k < num_tapers; k++) { for (int i = 0; i < n; i++) { tapered_signal[i] = signal->data[i] * sin(M_PI * (i + 0.5) / n); } fftw_execute(p); for (int i = 0; i < n / 2 + 1; i++) { psd[i] += (out[i][0] * out[i][0] + out[i][1] * out[i][1]) / n; } } for (int i = 0; i < n / 2 + 1; i++) { psd[i] /= num_tapers; } fftw_destroy_plan(p); fftw_free(tapered_signal); fftw_free(out); }