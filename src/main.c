/*
 * main.c: Main entry point for the EEG Signal Analysis Pipeline.
 *
 * This program reads a specified column from a CSV file, normalizes it,
 * and then performs a multi-stage comparative analysis:
 * 1. On the normalized original signal.
 * 2. On a pre-processed (detrended and filtered) version.
 * 3. On a uniformly quantized version of the filtered signal.
 * 4. On a non-uniformly (mu-law) quantized version of the filtered signal.
 * 5. A PCM binary encoding is generated from the uniform quantized signal for visualization.
 *
 * It benchmarks each analysis and writes all results for Python visualization.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "csv_parser.h"
#include "signal_processing.h"
#include "performance_monitor.h"

#define SAMPLING_RATE 256.0
#define QUANTIZATION_LEVELS 256 // 2^8 for 8-bit quantization
#define MU_PARAMETER 255.0      // Standard for 8-bit mu-law

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

    SignalData raw_signal = read_csv_column(filename, column_name);
    if (raw_signal.data == NULL) return 1;
    printf("Read %d data points from column '%s'.\n", raw_signal.count, column_name);

    printf("Normalizing signal to range [-1, 1]...\n");
    normalize_signal(&raw_signal);

    // --- Write config file for Python plotter ---
    int welch_windows[] = {128, 256, 512, 1024, 2048, 4096};
    int num_windows = sizeof(welch_windows) / sizeof(welch_windows[0]);
    snprintf(filepath, sizeof(filepath), "%s/config.txt", data_dir);
    FILE* config_file = fopen(filepath, "w");
    if (config_file) {
        fprintf(config_file, "sampling_rate,%.1f\n", SAMPLING_RATE);
        fprintf(config_file, "signal_length,%d\n", raw_signal.count);
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

    // --- Run 1: Analysis on NORMALIZED Original Signal ---
    printf("\n--- Running Analysis on NORMALIZED Signal ---\n");
    run_full_analysis(&raw_signal, "", data_dir, perf_file);

    // --- Create a copy for the filtered signal ---
    SignalData filtered_signal;
    filtered_signal.count = raw_signal.count;
    filtered_signal.data = (double*)malloc(raw_signal.count * sizeof(double));
    if (filtered_signal.data == NULL) { fprintf(stderr, "Failed to allocate memory for filtered signal.\n"); return 1; }
    memcpy(filtered_signal.data, raw_signal.data, raw_signal.count * sizeof(double));
    
    // --- Pre-processing Step 2: Detrend + Filter ---
    printf("\n--- Pre-processing Signal (Detrend + Filter) ---\n");
    detrend_signal(&filtered_signal);
    butterworth_bandpass_filter_4th_order(&filtered_signal, SAMPLING_RATE, 1.0, 40.0);

    // --- Run 2: Analysis on FILTERED Signal ---
    printf("\n--- Running Analysis on FILTERED Signal ---\n");
    run_full_analysis(&filtered_signal, "_filtered", data_dir, perf_file);

    // --- Create a copy for the uniformly quantized signal ---
    SignalData uniform_quant_signal;
    uniform_quant_signal.count = filtered_signal.count;
    uniform_quant_signal.data = (double*)malloc(filtered_signal.count * sizeof(double));
    if (uniform_quant_signal.data == NULL) { fprintf(stderr, "Failed to allocate memory for uniform quantized signal.\n"); return 1; }
    memcpy(uniform_quant_signal.data, filtered_signal.data, filtered_signal.count * sizeof(double));

    // --- Stage 3: Uniform Quantization ---
    printf("\n--- Uniformly Quantizing Filtered Signal (%d levels) ---\n", QUANTIZATION_LEVELS);
    quantize_signal_uniform(&uniform_quant_signal, QUANTIZATION_LEVELS);

    // --- Run 3: Analysis on UNIFORMLY QUANTIZED Signal ---
    printf("\n--- Running Analysis on UNIFORM QUANTIZED Signal ---\n");
    run_full_analysis(&uniform_quant_signal, "_uniform_quantized", data_dir, perf_file);

    // --- Create a copy for the mu-law quantized signal ---
    SignalData mu_law_quant_signal;
    mu_law_quant_signal.count = filtered_signal.count;
    mu_law_quant_signal.data = (double*)malloc(filtered_signal.count * sizeof(double));
    if (mu_law_quant_signal.data == NULL) { fprintf(stderr, "Failed to allocate memory for mu-law quantized signal.\n"); return 1; }
    memcpy(mu_law_quant_signal.data, filtered_signal.data, filtered_signal.count * sizeof(double));

    // --- Stage 4: Mu-Law Quantization ---
    printf("\n--- Mu-Law Quantizing Filtered Signal (mu=%.1f) ---\n", MU_PARAMETER);
    quantize_signal_mu_law(&mu_law_quant_signal, QUANTIZATION_LEVELS, MU_PARAMETER);

    // --- Run 4: Analysis on MU-LAW QUANTIZED Signal ---
    printf("\n--- Running Analysis on MU-LAW QUANTIZED Signal ---\n");
    run_full_analysis(&mu_law_quant_signal, "_mu_law_quantized", data_dir, perf_file);

    // --- Stage 5: PCM Encoding (for visualization only) ---
    printf("\n--- Generating PCM Encoded Stream for Visualization ---\n");
    char **pcm_encoded_stream = malloc(uniform_quant_signal.count * sizeof(char *));
    if (pcm_encoded_stream) {
        for (int i = 0; i < uniform_quant_signal.count; i++) {
            pcm_encoded_stream[i] = malloc(9 * sizeof(char)); // 8 bits + null terminator
        }
        pcm_encode(&uniform_quant_signal, pcm_encoded_stream, QUANTIZATION_LEVELS);
        
        snprintf(filepath, sizeof(filepath), "%s/pcm_encoded_stream.txt", data_dir);
        FILE *pcm_file = fopen(filepath, "w");
        if (pcm_file) {
            for (int i = 0; i < uniform_quant_signal.count; i++) {
                fprintf(pcm_file, "%s\n", pcm_encoded_stream[i]);
            }
            fclose(pcm_file);
        }
        for (int i = 0; i < uniform_quant_signal.count; i++) {
            free(pcm_encoded_stream[i]);
        }
        free(pcm_encoded_stream);
    }

    // --- Finalization ---
    fclose(perf_file);
    free_signal_data(&raw_signal);
    free_signal_data(&filtered_signal);
    free_signal_data(&uniform_quant_signal);
    free_signal_data(&mu_law_quant_signal);
    fftw_cleanup();

    printf("\nAnalysis complete for all stages.\n");
    return 0;
}

void run_full_analysis(SignalData* signal, const char* suffix, const char* data_dir, FILE* perf_file) {
    char filepath[256];
    char method_name[128];
    PerformanceMetrics metrics;
    struct rusage usage_start, usage_end;
    struct timespec timer_start;
    snprintf(filepath, sizeof(filepath), "%s/signal%s.txt", data_dir, suffix);
    write_array_to_file(filepath, signal->data, signal->count);
    printf("Calculating autocorrelation...\n");
    double *acf = (double *)malloc(signal->count * sizeof(double));
    if (!acf) { fprintf(stderr, "ACF allocation failed.\n"); return; }
    calculate_autocorrelation(signal, acf);
    snprintf(filepath, sizeof(filepath), "%s/acf%s.txt", data_dir, suffix);
    write_array_to_file(filepath, acf, signal->count);
    free(acf);
    printf("Analyzing PSD with Periodogram...\n");
    double* psd_periodogram_data = (double *)malloc((signal->count / 2 + 1) * sizeof(double));
    if (!psd_periodogram_data) { fprintf(stderr, "Periodogram allocation failed.\n"); return; }
    fftw_plan p = NULL;
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
    printf("Analyzing PSD with Welch (multiple windows)...\n");
    int welch_windows[] = {128, 256, 512, 1024, 2048, 4096};
    int num_windows = sizeof(welch_windows) / sizeof(welch_windows[0]);
    for (int i = 0; i < num_windows; i++) {
        int nperseg = welch_windows[i];
        if (signal->count < nperseg) continue;
        double* psd_welch_data = (double *)malloc((nperseg / 2 + 1) * sizeof(double));
        if (!psd_welch_data) { fprintf(stderr, "Welch allocation failed for window %d.\n", nperseg); continue; }
        start_timer(&timer_start); get_cpu_usage(&usage_start);
        psd_welch(signal, psd_welch_data, nperseg, nperseg / 2);
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