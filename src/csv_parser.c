#include "csv_parser.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LINE_LEN 2048
#define INITIAL_CAPACITY 50000

// Reads the header of the CSV to find the index of the target column
int get_column_index(const char *header, const char *column_name) {
    char *header_copy = strdup(header);
    char *token = strtok(header_copy, ",\n\r");
    int index = 0;
    while (token != NULL) {
        if (strcmp(token, column_name) == 0) {
            free(header_copy);
            return index;
        }
        token = strtok(NULL, ",\n\r");
        index++;
    }
    free(header_copy);
    return -1; // Not found
}

SignalData read_csv_column(const char *filename, const char *column_name) {
    SignalData result = {NULL, 0};
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        return result;
    }

    char line[MAX_LINE_LEN];

    if (fgets(line, sizeof(line), file) == NULL) {
        fprintf(stderr, "Error: Could not read header from CSV.\n");
        fclose(file);
        return result;
    }

    int col_index = get_column_index(line, column_name);
    if (col_index == -1) {
        fprintf(stderr, "Error: Column '%s' not found in '%s'.\n", column_name, filename);
        fclose(file);
        return result;
    }

    int capacity = INITIAL_CAPACITY;
    result.data = (double *)malloc(capacity * sizeof(double));
    if (!result.data) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        fclose(file);
        return result;
    }

    while (fgets(line, sizeof(line), file)) {
        if (result.count >= capacity) {
            capacity *= 2;
            double *temp = (double *)realloc(result.data, capacity * sizeof(double));
            if (!temp) {
                fprintf(stderr, "Error: Memory reallocation failed.\n");
                free(result.data);
                result.data = NULL;
                result.count = 0;
                fclose(file);
                return result;
            }
            result.data = temp;
        }

        char *token = strtok(line, ",\n\r");
        for (int i = 0; i < col_index; ++i) {
            token = strtok(NULL, ",\n\r");
        }

        if (token != NULL) {
            result.data[result.count++] = atof(token);
        }
    }

    fclose(file);
    return result;
}

void free_signal_data(SignalData *signal) {
    if (signal && signal->data) {
        free(signal->data);
        signal->data = NULL;
        signal->count = 0;
    }
}