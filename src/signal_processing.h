#ifndef SIGNAL_PROCESSING_H
#define SIGNAL_PROCESSING_H

#include "csv_parser.h"
#include <fftw3.h>

// --- Signal Generation Functions ---
void generate_gaussian_white_noise(SignalData *signal, int length, double amplitude);

// --- Pre-processing Functions ---
void normalize_signal(SignalData *signal);
void detrend_signal(SignalData *signal);
void butterworth_bandpass_filter_4th_order(SignalData *signal, double sample_rate, double low_cutoff, double high_cutoff);
void quantize_signal_uniform(SignalData *signal, int num_levels);
void quantize_signal_mu_law(SignalData *signal, int num_levels, double mu);
void pcm_encode(const SignalData *quantized_signal, char **encoded_output, int num_levels);


// --- Analysis Functions ---
// Original simple biased estimator
void calculate_autocorrelation_direct(const SignalData *signal, double *acf_result);
// NEW: Unbiased FFT-based estimator from reference code
void calculate_autocorrelation_fft(const SignalData *signal, double *acf_result);

// MODIFIED: Periodogram now uses a window option
void psd_periodogram(const SignalData *signal, double *psd, fftw_plan *p, int use_window);
// MODIFIED: Welch now performs local detrending and has robust scaling
void psd_welch_advanced(const SignalData *signal, double *psd, double fs, int nperseg, int noverlap);
// MODIFIED: Multitaper now uses Gaussian tapers
void psd_multitaper_gaussian(const SignalData *signal, double *psd, double fs, double nw, int k_tapers);

// --- Utility ---
void write_array_to_file(const char *filename, const double *data, int size);

#endif // SIGNAL_PROCESSING_H