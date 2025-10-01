#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "csv_parser.h"
#include "signal_processing.h"
#include "performance_monitor.h"

#define SAMPLING_RATE 256.0

// Global variables
static SignalData g_signal;
static double *g_psd_periodogram, *g_psd_welch, *g_psd_multitaper;
static fftw_plan g_periodogram_plan = NULL;
static int g_welch_nperseg;

// Wrappers
void run_periodogram_wrapper() { psd_periodogram(&g_signal, g_psd_periodogram, &g_periodogram_plan); }
void run_welch_wrapper() { psd_welch(&g_signal, g_psd_welch, g_welch_nperseg, g_welch_nperseg / 2); }
void run_multitaper_wrapper() { psd_multitaper(&g_signal, g_psd_multitaper); }

int main(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <csv_file> <column_name> <output_data_dir>\n", argv[0]);
        return 1;
    }

    const char *filename = argv[1];
    const char *column_name = argv[2];
    const char *data_dir = argv[3];
    char filepath[256];

    g_signal = read_csv_column(filename, column_name);
    if (g_signal.data == NULL) return 1;
    printf("Read %d data points.\n", g_signal.count);

    // --- Write config file for Python ---
    int welch_windows[] = {128, 256, 512, 1024, 2048, 4096};
    int num_windows = sizeof(welch_windows) / sizeof(welch_windows[0]);
    snprintf(filepath, sizeof(filepath), "%s/config.txt", data_dir);
    FILE* config_file = fopen(filepath, "w");
    if (config_file) {
        fprintf(config_file, "sampling_rate,%.1f\n", SAMPLING_RATE);
        fprintf(config_file, "signal_length,%d\n", g_signal.count);
        fprintf(config_file, "welch_windows");
        for (int i = 0; i < num_windows; i++) {
            fprintf(config_file, ",%d", welch_windows[i]);
        }
        fprintf(config_file, "\n");
        fclose(config_file);
    }

    // --- Autocorrelation ---
    printf("Calculating autocorrelation...\n");
    double *acf = (double *)malloc(g_signal.count * sizeof(double));
    calculate_autocorrelation(&g_signal, acf);
    snprintf(filepath, sizeof(filepath), "%s/acf.txt", data_dir);
    write_array_to_file(filepath, acf, g_signal.count);
    free(acf);
    
    // --- Performance Analysis ---
    FILE *perf_file = fopen("results/performance/performance.txt", "w");
    fprintf(perf_file, "Method,ExecutionTime_s,CPUUserTime_us,CPUSystemTime_us,PeakMemory_kb\n");
    PerformanceMetrics metrics;
    struct rusage usage_start, usage_end;
    struct timespec timer_start;

    // 1. Periodogram
    printf("Analyzing PSD with Periodogram...\n");
    g_psd_periodogram = (double *)malloc((g_signal.count / 2 + 1) * sizeof(double));
    start_timer(&timer_start); get_cpu_usage(&usage_start); run_periodogram_wrapper(); get_cpu_usage(&usage_end);
    metrics.execution_time_sec = stop_timer(&timer_start); metrics.peak_memory_kb = get_peak_memory_kb();
    metrics.cpu_user_time_us = (usage_end.ru_utime.tv_sec - usage_start.ru_utime.tv_sec) * 1000000L + (usage_end.ru_utime.tv_usec - usage_start.ru_utime.tv_usec);
    metrics.cpu_system_time_us = (usage_end.ru_stime.tv_sec - usage_start.ru_stime.tv_sec) * 1000000L + (usage_end.ru_stime.tv_usec - usage_start.ru_stime.tv_usec);
    fprintf(perf_file, "Periodogram,%.6f,%ld,%ld,%ld\n", metrics.execution_time_sec, metrics.cpu_user_time_us, metrics.cpu_system_time_us, metrics.peak_memory_kb);
    snprintf(filepath, sizeof(filepath), "%s/periodogram.txt", data_dir);
    write_array_to_file(filepath, g_psd_periodogram, g_signal.count / 2 + 1);
    free(g_psd_periodogram);

    // 2. Welch (Multiple Windows)
    printf("Analyzing PSD with Welch (multiple windows)...\n");
    for (int i = 0; i < num_windows; i++) {
        g_welch_nperseg = welch_windows[i];
        if (g_signal.count < g_welch_nperseg) {
            printf("  Skipping Welch with %d-point window (signal is too short).\n", g_welch_nperseg);
            continue;
        }
        printf("  Running Welch with %d-point window...\n", g_welch_nperseg);
        g_psd_welch = (double *)malloc((g_welch_nperseg / 2 + 1) * sizeof(double));
        start_timer(&timer_start); get_cpu_usage(&usage_start); run_welch_wrapper(); get_cpu_usage(&usage_end);
        metrics.execution_time_sec = stop_timer(&timer_start); metrics.peak_memory_kb = get_peak_memory_kb();
        metrics.cpu_user_time_us = (usage_end.ru_utime.tv_sec - usage_start.ru_utime.tv_sec) * 1000000L + (usage_end.ru_utime.tv_usec - usage_start.ru_utime.tv_usec);
        metrics.cpu_system_time_us = (usage_end.ru_stime.tv_sec - usage_start.ru_stime.tv_sec) * 1000000L + (usage_end.ru_stime.tv_usec - usage_start.ru_stime.tv_usec);
        fprintf(perf_file, "Welch_%d,%.6f,%ld,%ld,%ld\n", g_welch_nperseg, metrics.execution_time_sec, metrics.cpu_user_time_us, metrics.cpu_system_time_us, metrics.peak_memory_kb);
        snprintf(filepath, sizeof(filepath), "%s/welch_%d.txt", data_dir, g_welch_nperseg);
        write_array_to_file(filepath, g_psd_welch, g_welch_nperseg / 2 + 1);
        free(g_psd_welch);
    }

    // 3. Multitaper
    printf("Analyzing PSD with Multitaper...\n");
    g_psd_multitaper = (double *)malloc((g_signal.count / 2 + 1) * sizeof(double));
    start_timer(&timer_start); get_cpu_usage(&usage_start); run_multitaper_wrapper(); get_cpu_usage(&usage_end);
    metrics.execution_time_sec = stop_timer(&timer_start); metrics.peak_memory_kb = get_peak_memory_kb();
    metrics.cpu_user_time_us = (usage_end.ru_utime.tv_sec - usage_start.ru_utime.tv_sec) * 1000000L + (usage_end.ru_utime.tv_usec - usage_start.ru_utime.tv_usec);
    metrics.cpu_system_time_us = (usage_end.ru_stime.tv_sec - usage_start.ru_stime.tv_sec) * 1000000L + (usage_end.ru_stime.tv_usec - usage_start.ru_stime.tv_usec);
    fprintf(perf_file, "Multitaper,%.6f,%ld,%ld,%ld\n", metrics.execution_time_sec, metrics.cpu_user_time_us, metrics.cpu_system_time_us, metrics.peak_memory_kb);
    snprintf(filepath, sizeof(filepath), "%s/multitaper.txt", data_dir);
    write_array_to_file(filepath, g_psd_multitaper, g_signal.count / 2 + 1);
    free(g_psd_multitaper);

    // --- Finalization ---
    fclose(perf_file);
    snprintf(filepath, sizeof(filepath), "%s/signal.txt", data_dir);
    write_array_to_file(filepath, g_signal.data, g_signal.count);
    fftw_destroy_plan(g_periodogram_plan);
    free_signal_data(&g_signal);
    fftw_cleanup();

    printf("Analysis complete.\n");
    return 0;
}