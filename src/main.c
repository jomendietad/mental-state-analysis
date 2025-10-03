/*
 * main.c: Main entry point for the EEG Signal Analysis Pipeline.
 *
 * This program reads a specified column from a CSV file, performs a series of
 * signal processing analyses (Autocorrelation, Periodogram, Welch's, Multitaper)
 * on the raw signal, and then repeats the entire analysis on a pre-processed
 * version of the signal (detrended and band-pass filtered).
 *
 * It benchmarks the performance of each analysis method and writes all numerical
 * results to text files for later visualization by a Python script.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "csv_parser.h"
#include "signal_processing.h"
#include "performance_monitor.h"

#define SAMPLING_RATE 256.0

// A helper function to run the full analysis pipeline on a given signal
void run_full_analysis(SignalData* signal, const char* suffix, const char* data_dir, FILE* perf_file);

int main(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <csv_file> <column_name> <output_data_dir>\n", argv[0]);
        return 1;
    }

    const char *filename = argv[1];
    const char *column_name = argv[2];
    const char *data_dir = argv[3];
    char filepath[256];

    SignalData original_signal = read_csv_column(filename, column_name);
    if (original_signal.data == NULL) return 1;
    printf("Read %d data points from column '%s'.\n", original_signal.count, column_name);

    // --- Write config file for Python plotter ---
    int welch_windows[] = {128, 256, 512, 1024, 2048, 4096};
    int num_windows = sizeof(welch_windows) / sizeof(welch_windows[0]);
    snprintf(filepath, sizeof(filepath), "%s/config.txt", data_dir);
    FILE* config_file = fopen(filepath, "w");
    if (config_file) {
        fprintf(config_file, "sampling_rate,%.1f\n", SAMPLING_RATE);
        fprintf(config_file, "signal_length,%d\n", original_signal.count);
        fprintf(config_file, "welch_windows");
        for (int i = 0; i < num_windows; i++) { fprintf(config_file, ",%d", welch_windows[i]); }
        fprintf(config_file, "\n");
        fclose(config_file);
    }

    // --- Prepare performance results file ---
    snprintf(filepath, sizeof(filepath), "results/performance/performance.txt");
    FILE *perf_file = fopen(filepath, "w");
    if (!perf_file) {
        perror("Failed to open performance file");
        return 1;
    }
    fprintf(perf_file, "Method,ExecutionTime_s,CPUUserTime_us,CPUSystemTime_us,PeakMemory_kb\n");

    // --- Run 1: Analysis on Original Signal ---
    printf("\n--- Running Analysis on ORIGINAL Signal ---\n");
    run_full_analysis(&original_signal, "", data_dir, perf_file);

    // --- Create a copy of the signal to be filtered ---
    SignalData filtered_signal;
    filtered_signal.count = original_signal.count;
    filtered_signal.data = (double*)malloc(original_signal.count * sizeof(double));
    if (filtered_signal.data == NULL) {
        fprintf(stderr, "Failed to allocate memory for filtered signal.\n");
        return 1;
    }
    memcpy(filtered_signal.data, original_signal.data, original_signal.count * sizeof(double));
    
    // --- Pre-process the copied signal ---
    printf("\n--- Pre-processing Signal for Second Run (1-40 Hz 4th-Order Butterworth Band-pass) ---\n");
    detrend_signal(&filtered_signal);
    // MODIFIED: Using the new, more powerful 4th-order filter
    butterworth_bandpass_filter_4th_order(&filtered_signal, SAMPLING_RATE, 1.0, 40.0);

    // --- Run 2: Analysis on Filtered Signal ---
    printf("\n--- Running Analysis on FILTERED Signal ---\n");
    run_full_analysis(&filtered_signal, "_filtered", data_dir, perf_file);

    // --- Finalization ---
    fclose(perf_file);
    free_signal_data(&original_signal);
    free_signal_data(&filtered_signal);
    fftw_cleanup(); // Free resources used by FFTW

    printf("\nAnalysis complete for both original and filtered signals.\n");
    return 0;
}

/**
 * @brief Runs the complete set of analyses on a signal and records performance.
 *
 * This function calculates and saves the autocorrelation and three types of PSD
 * estimates (Periodogram, Welch, Multitaper). It measures the performance of
 * each PSD method and writes the results to the provided performance file.
 *
 * @param signal The input signal data.
 * @param suffix A string to append to output filenames (e.g., "" or "_filtered").
 * @param data_dir The directory to save the output data files.
 * @param perf_file A file pointer to the performance log.
 */
void run_full_analysis(SignalData* signal, const char* suffix, const char* data_dir, FILE* perf_file) {
    char filepath[256];
    char method_name[128];
    PerformanceMetrics metrics;
    struct rusage usage_start, usage_end;
    struct timespec timer_start;

    // Save the raw (or filtered) signal data
    snprintf(filepath, sizeof(filepath), "%s/signal%s.txt", data_dir, suffix);
    write_array_to_file(filepath, signal->data, signal->count);

    // --- Autocorrelation ---
    printf("Calculating autocorrelation...\n");
    double *acf = (double *)malloc(signal->count * sizeof(double));
    if (!acf) { fprintf(stderr, "ACF allocation failed.\n"); return; }
    calculate_autocorrelation(signal, acf);
    snprintf(filepath, sizeof(filepath), "%s/acf%s.txt", data_dir, suffix);
    write_array_to_file(filepath, acf, signal->count);
    free(acf);

    // --- PSD: Periodogram ---
    printf("Analyzing PSD with Periodogram...\n");
    double* psd_periodogram_data = (double *)malloc((signal->count / 2 + 1) * sizeof(double));
    if (!psd_periodogram_data) { fprintf(stderr, "Periodogram allocation failed.\n"); return; }
    fftw_plan p = NULL; // Plan is created inside the function
    start_timer(&timer_start); get_cpu_usage(&usage_start);
    psd_periodogram(signal, psd_periodogram_data, &p);
    get_cpu_usage(&usage_end); metrics.execution_time_sec = stop_timer(&timer_start);
    metrics.peak_memory_kb = get_peak_memory_kb();
    metrics.cpu_user_time_us = (usage_end.ru_utime.tv_sec - usage_start.ru_utime.tv_sec) * 1000000L + (usage_end.ru_utime.tv_usec - usage_start.ru_utime.tv_usec);
    metrics.cpu_system_time_us = (usage_end.ru_stime.tv_sec - usage_start.ru_stime.tv_sec) * 1000000L + (usage_end.ru_stime.tv_usec - usage_start.ru_stime.tv_usec);
    snprintf(method_name, sizeof(method_name), "Periodogram%s", suffix);
    fprintf(perf_file, "%s,%.6f,%ld,%ld,%ld\n", method_name, metrics.execution_time_sec, metrics.cpu_user_time_us, metrics.cpu_system_time_us, metrics.peak_memory_kb);
    snprintf(filepath, sizeof(filepath), "%s/periodogram%s.txt", data_dir, suffix);
    write_array_to_file(filepath, psd_periodogram_data, signal->count / 2 + 1);
    free(psd_periodogram_data);
    fftw_destroy_plan(p);

    // --- PSD: Welch's Method ---
    printf("Analyzing PSD with Welch (multiple windows)...\n");
    int welch_windows[] = {128, 256, 512, 1024, 2048, 4096};
    int num_windows = sizeof(welch_windows) / sizeof(welch_windows[0]);
    for (int i = 0; i < num_windows; i++) {
        int nperseg = welch_windows[i];
        if (signal->count < nperseg) continue;
        
        double* psd_welch_data = (double *)malloc((nperseg / 2 + 1) * sizeof(double));
        if (!psd_welch_data) { fprintf(stderr, "Welch allocation failed for window %d.\n", nperseg); continue; }
        
        start_timer(&timer_start); get_cpu_usage(&usage_start);
        psd_welch(signal, psd_welch_data, nperseg, nperseg / 2); // 50% overlap
        get_cpu_usage(&usage_end); metrics.execution_time_sec = stop_timer(&timer_start);
        metrics.peak_memory_kb = get_peak_memory_kb();
        metrics.cpu_user_time_us = (usage_end.ru_utime.tv_sec - usage_start.ru_utime.tv_sec) * 1000000L + (usage_end.ru_utime.tv_usec - usage_start.ru_utime.tv_usec);
        metrics.cpu_system_time_us = (usage_end.ru_stime.tv_sec - usage_start.ru_stime.tv_sec) * 1000000L + (usage_end.ru_stime.tv_usec - usage_start.ru_stime.tv_usec);
        
        snprintf(method_name, sizeof(method_name), "Welch_%d%s", nperseg, suffix);
        fprintf(perf_file, "%s,%.6f,%ld,%ld,%ld\n", method_name, metrics.execution_time_sec, metrics.cpu_user_time_us, metrics.cpu_system_time_us, metrics.peak_memory_kb);
        snprintf(filepath, sizeof(filepath), "%s/welch_%d%s.txt", data_dir, nperseg, suffix);
        write_array_to_file(filepath, psd_welch_data, nperseg / 2 + 1);
        free(psd_welch_data);
    }

    // --- PSD: Multitaper ---
    printf("Analyzing PSD with Multitaper...\n");
    double* psd_multitaper_data = (double *)malloc((signal->count / 2 + 1) * sizeof(double));
    if (!psd_multitaper_data) { fprintf(stderr, "Multitaper allocation failed.\n"); return; }
    start_timer(&timer_start); get_cpu_usage(&usage_start);
    psd_multitaper(signal, psd_multitaper_data);
    get_cpu_usage(&usage_end); metrics.execution_time_sec = stop_timer(&timer_start);
    metrics.peak_memory_kb = get_peak_memory_kb();
    metrics.cpu_user_time_us = (usage_end.ru_utime.tv_sec - usage_start.ru_utime.tv_sec) * 1000000L + (usage_end.ru_utime.tv_usec - usage_start.ru_utime.tv_usec);
    metrics.cpu_system_time_us = (usage_end.ru_stime.tv_sec - usage_start.ru_stime.tv_sec) * 1000000L + (usage_end.ru_stime.tv_usec - usage_start.ru_stime.tv_usec);
    snprintf(method_name, sizeof(method_name), "Multitaper%s", suffix);
    fprintf(perf_file, "%s,%.6f,%ld,%ld,%ld\n", method_name, metrics.execution_time_sec, metrics.cpu_user_time_us, metrics.cpu_system_time_us, metrics.peak_memory_kb);
    snprintf(filepath, sizeof(filepath), "%s/multitaper%s.txt", data_dir, suffix);
    write_array_to_file(filepath, psd_multitaper_data, signal->count / 2 + 1);
    free(psd_multitaper_data);
}