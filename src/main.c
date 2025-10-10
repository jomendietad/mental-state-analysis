#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "csv_parser.h"
#include "signal_processing.h"
#include "performance_monitor.h"

#define SAMPLING_RATE 256.0
#define MU_PARAMETER 255.0
#define WGN_AMPLITUDE 1.0
#define MULTITAPER_NW 2.5
#define MULTITAPER_K 4

void run_full_analysis(SignalData* signal, const char* suffix, const char* data_dir, FILE* perf_file);

int main(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <csv_file> <column_name> <output_data_dir>\n", argv[0]);
        return 1;
    }

    srand(time(NULL));

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
    
    printf("\n--- Pre-processing Signal (Detrend + Filter) ---\n");
    detrend_signal(&filtered_signal);
    butterworth_bandpass_filter_4th_order(&filtered_signal, SAMPLING_RATE, 1.0, 40.0);

    // --- Run 2: Analysis on FILTERED Signal ---
    printf("\n--- Running Analysis on FILTERED Signal ---\n");
    run_full_analysis(&filtered_signal, "_filtered", data_dir, perf_file);

    // --- Stage 3: Uniform Quantization (8-bit) ---
    SignalData uniform_quant_signal;
    uniform_quant_signal.count = filtered_signal.count;
    uniform_quant_signal.data = (double*)malloc(filtered_signal.count * sizeof(double));
    if (uniform_quant_signal.data == NULL) { fprintf(stderr, "Failed to allocate memory for uniform quantized signal.\n"); return 1; }
    memcpy(uniform_quant_signal.data, filtered_signal.data, filtered_signal.count * sizeof(double));
    printf("\n--- Uniformly Quantizing Filtered Signal (8-bit, 256 levels) ---\n");
    quantize_signal_uniform(&uniform_quant_signal, 256);
    printf("\n--- Running Analysis on UNIFORM QUANTIZED Signal ---\n");
    run_full_analysis(&uniform_quant_signal, "_uniform_quantized", data_dir, perf_file);

    // --- Stage 4: Mu-Law Quantization (8-bit) ---
    SignalData mu_law_quant_signal;
    mu_law_quant_signal.count = filtered_signal.count;
    mu_law_quant_signal.data = (double*)malloc(filtered_signal.count * sizeof(double));
    if (mu_law_quant_signal.data == NULL) { fprintf(stderr, "Failed to allocate memory for mu-law quantized signal.\n"); return 1; }
    memcpy(mu_law_quant_signal.data, filtered_signal.data, filtered_signal.count * sizeof(double));
    printf("\n--- Mu-Law Quantizing Filtered Signal (mu=%.1f) ---\n", MU_PARAMETER);
    quantize_signal_mu_law(&mu_law_quant_signal, 256, MU_PARAMETER);
    printf("\n--- Running Analysis on MU-LAW QUANTIZED Signal ---\n");
    run_full_analysis(&mu_law_quant_signal, "_mu_law_quantized", data_dir, perf_file);

    // --- Generate lower-bit-rate signals for error analysis only ---
    printf("\n--- Generating lower-bit-rate signals for error analysis ---\n");
    // 4-bit
    SignalData quant_4bit_signal;
    quant_4bit_signal.count = filtered_signal.count;
    quant_4bit_signal.data = (double*)malloc(filtered_signal.count * sizeof(double));
    if (quant_4bit_signal.data) {
        memcpy(quant_4bit_signal.data, filtered_signal.data, filtered_signal.count * sizeof(double));
        quantize_signal_uniform(&quant_4bit_signal, 16); // 2^4 = 16 levels
        snprintf(filepath, sizeof(filepath), "%s/signal_quantized_4bit.txt", data_dir);
        write_array_to_file(filepath, quant_4bit_signal.data, quant_4bit_signal.count);
        free_signal_data(&quant_4bit_signal);
    }
    // 2-bit
    SignalData quant_2bit_signal;
    quant_2bit_signal.count = filtered_signal.count;
    quant_2bit_signal.data = (double*)malloc(filtered_signal.count * sizeof(double));
    if (quant_2bit_signal.data) {
        memcpy(quant_2bit_signal.data, filtered_signal.data, filtered_signal.count * sizeof(double));
        quantize_signal_uniform(&quant_2bit_signal, 4); // 2^2 = 4 levels
        snprintf(filepath, sizeof(filepath), "%s/signal_quantized_2bit.txt", data_dir);
        write_array_to_file(filepath, quant_2bit_signal.data, quant_2bit_signal.count);
        free_signal_data(&quant_2bit_signal);
    }

    // --- PCM Encoding (from 8-bit signal for visualization) ---
    printf("\n--- Generating PCM Encoded Stream for Visualization ---\n");
    char **pcm_encoded_stream = malloc(uniform_quant_signal.count * sizeof(char *));
    if (pcm_encoded_stream) {
        for (int i = 0; i < uniform_quant_signal.count; i++) {
            pcm_encoded_stream[i] = malloc(9 * sizeof(char));
        }
        pcm_encode(&uniform_quant_signal, pcm_encoded_stream, 256);
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

    // --- Autocorrelation (Unbiased FFT method) ---
    printf("Calculating autocorrelation (FFT method)...\n");
    double *acf_fft = (double *)malloc(signal->count * sizeof(double));
    if (acf_fft) {
        calculate_autocorrelation_fft(signal, acf_fft);
        snprintf(filepath, sizeof(filepath), "%s/acf_fft%s.txt", data_dir, suffix);
        write_array_to_file(filepath, acf_fft, signal->count);
        free(acf_fft);
    }

    // --- PSD: Periodogram ---
    printf("Analyzing PSD with Periodogram...\n");
    double* psd_periodogram_data = (double *)malloc((signal->count / 2 + 1) * sizeof(double));
    if (!psd_periodogram_data) { fprintf(stderr, "Periodogram allocation failed.\n"); return; }
    fftw_plan p = NULL;
    start_timer(&timer_start); get_cpu_usage(&usage_start);
    psd_periodogram(signal, psd_periodogram_data, &p, 1); // Use window
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

    // --- PSD: Welch's Method (Advanced) ---
    printf("Analyzing PSD with Welch (multiple windows)...\n");
    int welch_windows[] = {128, 256, 512, 1024, 2048, 4096};
    int num_windows = sizeof(welch_windows) / sizeof(welch_windows[0]);
    for (int i = 0; i < num_windows; i++) {
        int nperseg = welch_windows[i];
        if (signal->count < nperseg) continue;
        double* psd_welch_data = (double *)malloc((nperseg / 2 + 1) * sizeof(double));
        if (!psd_welch_data) { fprintf(stderr, "Welch allocation failed for window %d.\n", nperseg); continue; }
        start_timer(&timer_start); get_cpu_usage(&usage_start);
        psd_welch_advanced(signal, psd_welch_data, SAMPLING_RATE, nperseg, nperseg / 2);
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

    // --- PSD: Multitaper (Gaussian) ---
    printf("Analyzing PSD with Multitaper (Gaussian)...\n");
    double* psd_multitaper_data = (double *)malloc((signal->count / 2 + 1) * sizeof(double));
    if (!psd_multitaper_data) { fprintf(stderr, "Multitaper allocation failed.\n"); return; }
    start_timer(&timer_start); get_cpu_usage(&usage_start);
    psd_multitaper_gaussian(signal, psd_multitaper_data, SAMPLING_RATE, MULTITAPER_NW, MULTITAPER_K);
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