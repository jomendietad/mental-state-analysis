#ifndef CSV_PARSER_H
#define CSV_PARSER_H

// Structure to hold the parsed data
typedef struct {
    double *data;
    int count;
} SignalData;

// Function to read a specific column from a CSV file
SignalData read_csv_column(const char *filename, const char *column_name);

// Function to free the memory allocated for SignalData
void free_signal_data(SignalData *signal);

#endif // CSV_PARSER_H