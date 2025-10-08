#include "signal_processing.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>

// --- Pre-processing Implementations ---

void normalize_signal(SignalData *signal) {
    // ... (code is unchanged)
    if (signal->count == 0) return;
    double max_abs = 0.0;
    for (int i = 0; i < signal->count; i++) {
        if (fabs(signal->data[i]) > max_abs) max_abs = fabs(signal->data[i]);
    }
    if (max_abs == 0.0) return;
    for (int i = 0; i < signal->count; i++) {
        signal->data[i] /= max_abs;
    }
}

void detrend_signal(SignalData *signal) {
    // ... (code is unchanged)
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

void butterworth_bandpass_filter_4th_order(SignalData *signal, double sample_rate, double low_cutoff, double high_cutoff) {
    // ... (code is unchanged)
    int n = signal->count;
    if (n < 5) return;
    const double b[] = {0.0469766, 0, -0.0939532, 0, 0.0469766};
    const double a[] = {1.0, -3.0348, 3.8155, -2.2855, 0.5455};
    double *temp_data = (double *)malloc(n * sizeof(double));
    if (!temp_data) { fprintf(stderr, "Memory allocation failed for filter buffer.\n"); return; }
    memcpy(temp_data, signal->data, n * sizeof(double));
    double x[5] = {0}; double y[5] = {0};
    for (int i = 0; i < n; i++) {
        memmove(&x[1], &x[0], 4 * sizeof(double));
        memmove(&y[1], &y[0], 4 * sizeof(double));
        x[0] = temp_data[i];
        y[0] = b[0] * x[0] + b[1] * x[1] + b[2] * x[2] + b[3] * x[3] + b[4] * x[4]
                         - a[1] * y[1] - a[2] * y[2] - a[3] * y[3] - a[4] * y[4];
        signal->data[i] = y[0];
    }
    free(temp_data);
}

// RENAMED for clarity
void quantize_signal_uniform(SignalData *signal, int num_levels) {
    // ... (code is unchanged)
    if (signal->count == 0 || num_levels <= 1) return;
    double min_val = DBL_MAX;
    double max_val = -DBL_MAX;
    for (int i = 0; i < signal->count; i++) {
        if (signal->data[i] < min_val) min_val = signal->data[i];
        if (signal->data[i] > max_val) max_val = signal->data[i];
    }
    if (max_val == min_val) return;
    double delta = (max_val - min_val) / (double)(num_levels - 1);
    for (int i = 0; i < signal->count; i++) {
        double level_index = round((signal->data[i] - min_val) / delta);
        signal->data[i] = min_val + level_index * delta;
    }
}

// --- NEW: Non-Uniform Quantization using mu-Law Algorithm ---
// Helper function for mu-law compression
double mu_law_compress(double x, double mu) {
    return copysign(1.0, x) * (log(1.0 + mu * fabs(x)) / log(1.0 + mu));
}

// Helper function for mu-law expansion
double mu_law_expand(double y, double mu) {
    return copysign(1.0, y) * (1.0 / mu) * (pow(1.0 + mu, fabs(y)) - 1.0);
}

void quantize_signal_mu_law(SignalData *signal, int num_levels, double mu) {
    if (signal->count == 0 || num_levels <= 1) {
        return;
    }
    // Note: This function assumes the input signal is already normalized to [-1, 1]

    // 1. Compress the signal using mu-law formula
    for (int i = 0; i < signal->count; i++) {
        signal->data[i] = mu_law_compress(signal->data[i], mu);
    }

    // 2. Apply UNIFORM quantization on the COMPRESSED signal
    // The compressed signal theoretically ranges from -1 to 1.
    double min_val = -1.0;
    double max_val = 1.0;
    double delta = (max_val - min_val) / (double)(num_levels - 1);

    for (int i = 0; i < signal->count; i++) {
        double level_index = round((signal->data[i] - min_val) / delta);
        signal->data[i] = min_val + level_index * delta;
    }

    // 3. Expand the quantized signal back to the original domain
    for (int i = 0; i < signal->count; i++) {
        signal->data[i] = mu_law_expand(signal->data[i], mu);
    }
}


// --- Analysis and Utility Functions (Unchanged) ---
// ... (The rest of the functions are identical)
void write_array_to_file(const char *filename, const double *data, int size) { FILE *f = fopen(filename, "w"); if (!f) { perror("Failed to open file"); return; } for (int i = 0; i < size; i++) { fprintf(f, "%.8f\n", data[i]); } fclose(f); }
void calculate_autocorrelation(const SignalData *signal, double *acf_result) { int n = signal->count; for (int lag = 0; lag < n; lag++) { double sum = 0.0; for (int i = 0; i < n - lag; i++) { sum += signal->data[i] * signal->data[i + lag]; } acf_result[lag] = sum / n; } }
void psd_periodogram(const SignalData *signal, double *psd, fftw_plan *p) { int n = signal->count; fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (n / 2 + 1)); if (*p == NULL) { *p = fftw_plan_dft_r2c_1d(n, (double*)signal->data, out, FFTW_ESTIMATE); } fftw_execute(*p); for (int i = 0; i < n / 2 + 1; i++) { psd[i] = (out[i][0] * out[i][0] + out[i][1] * out[i][1]) / n; } fftw_free(out); }
void psd_welch(const SignalData *signal, double *psd, int nperseg, int noverlap) { int n = signal->count; int step = nperseg - noverlap; if (step <= 0) return; int num_segments = (n - noverlap) / step; if (num_segments == 0) { return; } double *segment = (double*) fftw_malloc(sizeof(double) * nperseg); fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (nperseg / 2 + 1)); double *window = (double*) malloc(sizeof(double) * nperseg); for(int i = 0; i < nperseg; i++) { window[i] = 0.5 * (1 - cos(2 * M_PI * i / (nperseg - 1))); } fftw_plan p = fftw_plan_dft_r2c_1d(nperseg, segment, out, FFTW_ESTIMATE); memset(psd, 0, sizeof(double) * (nperseg / 2 + 1)); for (int i = 0; i < num_segments; i++) { for (int j = 0; j < nperseg; j++) { segment[j] = signal->data[i * step + j] * window[j]; } fftw_execute(p); for (int j = 0; j < nperseg / 2 + 1; j++) { psd[j] += (out[j][0] * out[j][0] + out[j][1] * out[j][1]) / nperseg; } } for (int i = 0; i < nperseg / 2 + 1; i++) { psd[i] /= num_segments; } fftw_destroy_plan(p); fftw_free(segment); fftw_free(out); free(window); }
void psd_multitaper(const SignalData *signal, double *psd) { int n = signal->count; int num_tapers = 4; double *tapered_signal = (double*) fftw_malloc(sizeof(double) * n); fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (n / 2 + 1)); fftw_plan p = fftw_plan_dft_r2c_1d(n, tapered_signal, out, FFTW_ESTIMATE); memset(psd, 0, sizeof(double) * (n / 2 + 1)); for (int k = 0; k < num_tapers; k++) { for (int i = 0; i < n; i++) { tapered_signal[i] = signal->data[i] * sin(M_PI * (i + 0.5) / n); } fftw_execute(p); for (int i = 0; i < n / 2 + 1; i++) { psd[i] += (out[i][0] * out[i][0] + out[i][1] * out[i][1]) / n; } } for (int i = 0; i < n / 2 + 1; i++) { psd[i] /= num_tapers; } fftw_destroy_plan(p); fftw_free(tapered_signal); fftw_free(out); }