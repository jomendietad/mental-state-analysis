#define _GNU_SOURCE
#include "performance_monitor.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void start_timer(struct timespec *start) {
    clock_gettime(CLOCK_MONOTONIC, start);
}

double stop_timer(struct timespec *start) {
    struct timespec end;
    clock_gettime(CLOCK_MONOTONIC, &end);
    return (end.tv_sec - start->tv_sec) + (end.tv_nsec - start->tv_nsec) / 1e9;
}

long get_peak_memory_kb() {
    FILE* fp = fopen("/proc/self/status", "r");
    if (!fp) {
        return -1;
    }
    long peak_mem = -1;
    char line[128];
    while (fgets(line, sizeof(line), fp)) {
        if (strncmp(line, "VmPeak:", 7) == 0) {
            sscanf(line, "VmPeak: %ld kB", &peak_mem);
            break;
        }
    }
    fclose(fp);
    return peak_mem;
}

void get_cpu_usage(struct rusage *usage) {
    if (getrusage(RUSAGE_SELF, usage) != 0) {
        perror("getrusage");
    }
}