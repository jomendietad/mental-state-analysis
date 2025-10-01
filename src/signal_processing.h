#ifndef SIGNAL_PROCESSING_H
#define SIGNAL_PROCESSING_H

#include "csv_parser.h"
#include <fftw3.h>

// --- Analysis Functions ---
void calculate_autocorrelation(const SignalData *signal, double *acf_result);
void psd_periodogram(const SignalData *signal, double *psd, fftw_plan *p);
void psd_welch(const SignalData *signal, double *psd, int nperseg, int noverlap);
void psd_multitaper(const SignalData *signal, double *psd);

// --- Utility ---
void write_array_to_file(const char *filename, const double *data, int size);

#endif // SIGNAL_PROCESSING_H