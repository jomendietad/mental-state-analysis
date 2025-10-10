#include "signal_processing.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <time.h>

// --- Signal Generation ---
void generate_gaussian_white_noise(SignalData *signal, int length, double amplitude) {
    if (length <= 0) return;

    signal->data = (double *)malloc(length * sizeof(double));
    if (!signal->data) {
        fprintf(stderr, "Memory allocation failed for GWN signal.\n");
        signal->count = 0;
        return;
    }
    signal->count = length;

    // The power of the noise will be amplitude^2
    for (int i = 0; i < length; i += 2) {
        double u1 = ((double)rand() + 1.0) / ((double)RAND_MAX + 2.0);
        double u2 = ((double)rand() + 1.0) / ((double)RAND_MAX + 2.0);

        double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
        double z1 = sqrt(-2.0 * log(u1)) * sin(2.0 * M_PI * u2);

        signal->data[i] = amplitude * z0;
        if (i + 1 < length) {
            signal->data[i + 1] = amplitude * z1;
        }
    }
}

// --- Pre-processing Implementations ---

void normalize_signal(SignalData *signal) {
    if (signal->count == 0) return;
    double max_abs = 0.0;
    for (int i = 0; i < signal->count; i++) {
        if (fabs(signal->data[i]) > max_abs) {
            max_abs = fabs(signal->data[i]);
        }
    }
    if (max_abs == 0.0) return;
    for (int i = 0; i < signal->count; i++) {
        signal->data[i] /= max_abs;
    }
}

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

void butterworth_bandpass_filter_4th_order(SignalData *signal, double sample_rate, double low_cutoff, double high_cutoff) {
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

void quantize_signal_uniform(SignalData *signal, int num_levels) {
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

double mu_law_compress(double x, double mu) {
    return copysign(1.0, x) * (log(1.0 + mu * fabs(x)) / log(1.0 + mu));
}
double mu_law_expand(double y, double mu) {
    return copysign(1.0, y) * (1.0 / mu) * (pow(1.0 + mu, fabs(y)) - 1.0);
}
void quantize_signal_mu_law(SignalData *signal, int num_levels, double mu) {
    if (signal->count == 0 || num_levels <= 1) return;
    for (int i = 0; i < signal->count; i++) {
        signal->data[i] = mu_law_compress(signal->data[i], mu);
    }
    double delta = 2.0 / (double)(num_levels - 1);
    for (int i = 0; i < signal->count; i++) {
        double level_index = round((signal->data[i] - (-1.0)) / delta);
        signal->data[i] = -1.0 + level_index * delta;
    }
    for (int i = 0; i < signal->count; i++) {
        signal->data[i] = mu_law_expand(signal->data[i], mu);
    }
}

void int_to_binary_string(int n, char *str) {
    str[8] = '\0';
    for (int i = 7; i >= 0; i--) {
        str[i] = (n & 1) + '0';
        n >>= 1;
    }
}

void pcm_encode(const SignalData *quantized_signal, char **encoded_output, int num_levels) {
    if (quantized_signal->count == 0 || num_levels <= 1) return;
    double min_val = DBL_MAX;
    double max_val = -DBL_MAX;
    for (int i = 0; i < quantized_signal->count; i++) {
        if (quantized_signal->data[i] < min_val) min_val = quantized_signal->data[i];
        if (quantized_signal->data[i] > max_val) max_val = quantized_signal->data[i];
    }
    if (max_val == min_val) {
        for(int i = 0; i < quantized_signal->count; i++) {
            int_to_binary_string(0, encoded_output[i]);
        }
        return;
    }
    double delta = (max_val - min_val) / (double)(num_levels - 1);
    for (int i = 0; i < quantized_signal->count; i++) {
        int level_index = (int)round((quantized_signal->data[i] - min_val) / delta);
        int_to_binary_string(level_index, encoded_output[i]);
    }
}

// --- Analysis Functions ---

// Original simple (biased) estimator
void calculate_autocorrelation_direct(const SignalData *signal, double *acf_result) {
    int n = signal->count;
    for (int lag = 0; lag < n; lag++) {
        double sum = 0.0;
        for (int i = 0; i < n - lag; i++) {
            sum += signal->data[i] * signal->data[i + lag];
        }
        acf_result[lag] = sum / n;
    }
}

// Unbiased FFT-based estimator (based on AutoCorr4 and unbiased_autocorrelation)
void calculate_autocorrelation_fft(const SignalData *signal, double *acf_result) {
    int n = signal->count;
    if (n == 0) return;

    // Demean the signal
    double mean = 0.0;
    for (int i = 0; i < n; i++) mean += signal->data[i];
    mean /= n;
    
    double *demeaned_signal = (double *)malloc(n * sizeof(double));
    if (!demeaned_signal) return;
    for (int i = 0; i < n; i++) demeaned_signal[i] = signal->data[i] - mean;

    // Zero pad to next power of two for speed (like AutoCorr4)
    int padded_size = 1;
    while (padded_size < 2 * n - 1) padded_size *= 2;

    double *padded_input = (double*) fftw_malloc(sizeof(double) * padded_size);
    fftw_complex *fft_output = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (padded_size / 2 + 1));
    double *ifft_output = (double*) fftw_malloc(sizeof(double) * padded_size);

    if (!padded_input || !fft_output || !ifft_output) {
        if (padded_input) fftw_free(padded_input);
        if (fft_output) fftw_free(fft_output);
        if (ifft_output) fftw_free(ifft_output);
        free(demeaned_signal);
        return;
    }

    memcpy(padded_input, demeaned_signal, n * sizeof(double));
    for (int i = n; i < padded_size; i++) padded_input[i] = 0.0;

    fftw_plan p_fwd = fftw_plan_dft_r2c_1d(padded_size, padded_input, fft_output, FFTW_ESTIMATE);
    fftw_plan p_bwd = fftw_plan_dft_c2r_1d(padded_size, fft_output, ifft_output, FFTW_ESTIMATE);

    fftw_execute(p_fwd);

    // Multiply by complex conjugate: |X(f)|^2
    for (int i = 0; i < (padded_size / 2 + 1); i++) {
        double real = fft_output[i][0];
        double imag = fft_output[i][1];
        fft_output[i][0] = real * real + imag * imag;
        fft_output[i][1] = 0.0;
    }

    fftw_execute(p_bwd);

    // Normalize for unbiased estimate
    for (int i = 0; i < n; i++) {
        if (n - i > 0) {
            // The IFFT result from FFTW is unnormalized, so we divide by padded_size
            acf_result[i] = (ifft_output[i] / padded_size) / (n - i);
        } else {
            acf_result[i] = 0;
        }
    }

    fftw_destroy_plan(p_fwd);
    fftw_destroy_plan(p_bwd);
    fftw_free(padded_input);
    fftw_free(fft_output);
    fftw_free(ifft_output);
    free(demeaned_signal);
}

// Periodogram with a window option
void psd_periodogram(const SignalData *signal, double *psd, fftw_plan *p, int use_window) {
    int n = signal->count;
    double *input_data = (double*)signal->data;
    double *windowed_data = NULL;
    double scale = n;

    if (use_window) {
        windowed_data = (double *)malloc(n * sizeof(double));
        if (!windowed_data) return;
        
        double window_sum_sq = 0.0;
        for (int i = 0; i < n; i++) {
            double window_val = 0.5 * (1 - cos(2 * M_PI * i / (n - 1))); // Hanning
            windowed_data[i] = signal->data[i] * window_val;
            window_sum_sq += window_val * window_val;
        }
        input_data = windowed_data;
        scale = window_sum_sq;
    }

    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (n / 2 + 1));
    if (*p == NULL) {
        *p = fftw_plan_dft_r2c_1d(n, input_data, out, FFTW_ESTIMATE);
    }
    // Use the specific execution function for the plan
    fftw_execute_dft_r2c(*p, input_data, out);

    for (int i = 0; i < n / 2 + 1; i++) {
        psd[i] = (out[i][0] * out[i][0] + out[i][1] * out[i][1]) / scale;
    }

    if (windowed_data) free(windowed_data);
    fftw_free(out);
}

// Welch with local detrending and robust scaling
void psd_welch_advanced(const SignalData *signal, double *psd, double fs, int nperseg, int noverlap) {
    int n = signal->count;
    int step = nperseg - noverlap;
    if (step <= 0 || n < nperseg) return;

    int num_segments = (n - noverlap) / step;
    if (num_segments == 0) return;

    double *segment = (double*) fftw_malloc(sizeof(double) * nperseg);
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (nperseg / 2 + 1));
    double *window = (double*) malloc(sizeof(double) * nperseg);
    double window_energy = 0.0;
    for(int i = 0; i < nperseg; i++) {
        window[i] = 0.5 * (1 - cos(2 * M_PI * i / (nperseg - 1))); // Hanning
        window_energy += window[i] * window[i];
    }

    fftw_plan p = fftw_plan_dft_r2c_1d(nperseg, segment, out, FFTW_ESTIMATE);
    memset(psd, 0, sizeof(double) * (nperseg / 2 + 1));

    for (int i = 0; i < num_segments; i++) {
        // Local detrending (mean removal)
        double mean = 0.0;
        for (int j = 0; j < nperseg; j++) mean += signal->data[i * step + j];
        mean /= nperseg;

        for (int j = 0; j < nperseg; j++) {
            segment[j] = (signal->data[i * step + j] - mean) * window[j];
        }
        
        fftw_execute(p);
        
        for (int j = 0; j < nperseg / 2 + 1; j++) {
            psd[j] += (out[j][0] * out[j][0] + out[j][1] * out[j][1]);
        }
    }

    // Averaging and scaling
    double scale = fs * window_energy;
    for (int i = 0; i < nperseg / 2 + 1; i++) {
        psd[i] /= num_segments;
        psd[i] /= scale;
        if (i > 0 && i < nperseg / 2) { // Apply one-sided scaling
            psd[i] *= 2.0;
        }
    }

    fftw_destroy_plan(p);
    fftw_free(segment);
    fftw_free(out);
    free(window);
}

// Multitaper with Gaussian tapers
void psd_multitaper_gaussian(const SignalData *signal, double *psd, double fs, double nw, int k_tapers) {
    int n = signal->count;
    if (n < 2 || k_tapers <= 0) return;

    double *tapered_signal = (double*) fftw_malloc(sizeof(double) * n);
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (n / 2 + 1));
    fftw_plan p = fftw_plan_dft_r2c_1d(n, tapered_signal, out, FFTW_ESTIMATE);
    
    memset(psd, 0, sizeof(double) * (n / 2 + 1));
    
    double sigma = n / (2.0 * M_PI * nw);

    for (int k = 0; k < k_tapers; k++) {
        double center_shift = (k - (k_tapers - 1.0) / 2.0) * n / (k_tapers * 5.0);
        double norm_factor = 0.0;

        // Generate and apply taper
        for (int i = 0; i < n; i++) {
            double n_centered = (i - n / 2.0 - center_shift) / sigma;
            double taper_val = exp(-0.5 * n_centered * n_centered);
            tapered_signal[i] = signal->data[i] * taper_val;
            norm_factor += taper_val * taper_val;
        }

        fftw_execute(p);

        // Accumulate PSD, normalizing by taper energy and fs
        double scale = fs * norm_factor;
        for (int i = 0; i < n / 2 + 1; i++) {
            double power = (out[i][0] * out[i][0] + out[i][1] * out[i][1]) / scale;
            if (i > 0 && i < n / 2) power *= 2.0;
            psd[i] += power;
        }
    }

    // Average over the number of tapers
    for (int i = 0; i < n / 2 + 1; i++) {
        psd[i] /= k_tapers;
    }

    fftw_destroy_plan(p);
    fftw_free(tapered_signal);
    fftw_free(out);
}

// --- Utility ---
void write_array_to_file(const char *filename, const double *data, int size) {
    FILE *f = fopen(filename, "w");
    if (!f) {
        perror("Failed to open file");
        return;
    }
    for (int i = 0; i < size; i++) {
        fprintf(f, "%.8f\n", data[i]);
    }
    fclose(f);
}