#ifndef SIGNAL_PROCESSING_H
#define SIGNAL_PROCESSING_H

#include "csv_parser.h"
#include <fftw3.h>

// --- Pre-processing Functions ---
void detrend_signal(SignalData *signal);
// NEW: A superior, 4th-order Butterworth filter function
void butterworth_bandpass_filter_4th_order(SignalData *signal, double sample_rate, double low_cutoff, double high_cutoff);


// --- Analysis Functions ---
void calculate_autocorrelation(const SignalData *signal, double *acf_result);
void psd_periodogram(const SignalData *signal, double *psd, fftw_plan *p);
void psd_welch(const SignalData *signal, double *psd, int nperseg, int noverlap);
void psd_multitaper(const SignalData *signal, double *psd);

// --- Utility ---
void write_array_to_file(const char *filename, const double *data, int size);

#endif // SIGNAL_PROCESSING_H