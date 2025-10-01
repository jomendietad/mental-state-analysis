#ifndef PERFORMANCE_MONITOR_H
#define PERFORMANCE_MONITOR_H

#include <time.h>
#include <sys/resource.h>

// Structure to hold performance metrics
typedef struct {
    double execution_time_sec;
    long cpu_user_time_us;
    long cpu_system_time_us;
    long peak_memory_kb;
} PerformanceMetrics;

// Function prototypes
void start_timer(struct timespec *start);
double stop_timer(struct timespec *start);
long get_peak_memory_kb();
void get_cpu_usage(struct rusage *usage);

#endif // PERFORMANCE_MONITOR_H