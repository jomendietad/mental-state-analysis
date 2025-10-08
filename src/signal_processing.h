#ifndef SIGNAL_PROCESSING_H
#define SIGNAL_PROCESSING_H

#include "csv_parser.h"
#include <fftw3.h>

// --- Pre-processing Functions ---
void normalize_signal(SignalData *signal);
void detrend_signal(SignalData *signal);
void butterworth_bandpass_filter_4th_order(SignalData *signal, double sample_rate, double low_cutoff, double high_cutoff);
void quantize_signal_uniform(SignalData *signal, int num_levels);
void quantize_signal_mu_law(SignalData *signal, int num_levels, double mu);
// NEW: PCM Encoding function
void pcm_encode(const SignalData *quantized_signal, char **encoded_output, int num_levels);


// --- Analysis Functions ---
void calculate_autocorrelation(const SignalData *signal, double *acf_result);
void psd_periodogram(const SignalData *signal, double *psd, fftw_plan *p);
void psd_welch(const SignalData *signal, double *psd, int nperseg, int noverlap);
void psd_multitaper(const SignalData *signal, double *psd);

// --- Utility ---
void write_array_to_file(const char *filename, const double *data, int size);

#endif // SIGNAL_PROCESSING_H