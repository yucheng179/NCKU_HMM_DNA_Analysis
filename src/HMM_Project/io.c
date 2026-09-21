#include "io.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

char *read_sequence_from_file(const char *filename, int *len) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    int capacity = 1024;
    char *sequence = malloc(capacity * sizeof(char));
    if (!sequence) {
        perror("Memory allocation failed");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    *len = 0;
    char c;
    while ((c = fgetc(file)) != EOF) {
        if (tolower(c) == 'a' || tolower(c) == 'c' || tolower(c) == 'g' || tolower(c) == 't') {
            if (*len >= capacity) {
                capacity *= 2;
                sequence = realloc(sequence, capacity * sizeof(char));
                if (!sequence) {
                    perror("Memory reallocation failed");
                    fclose(file);
                    exit(EXIT_FAILURE);
                }
            }
            sequence[(*len)++] = c;
        }
    }

    fclose(file);

    if (*len == 0) {
        fprintf(stderr, "Error: The file %s contains no valid DNA sequence.", filename);
        free(sequence);
        exit(EXIT_FAILURE);
    }

    sequence[*len] = '\0';
    return sequence;
}