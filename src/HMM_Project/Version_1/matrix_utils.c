#include "matrix_utils.h"
#include "hmm.h" 
#include <stdio.h>

void calculate_state_transition_matrix(char *sequence, int *hidden_states, int seq_len,
                                       double high_matrix[OBSERVATION_STATE][OBSERVATION_STATE],
                                       double low_matrix[OBSERVATION_STATE][OBSERVATION_STATE]) {
    int high_counts[OBSERVATION_STATE][OBSERVATION_STATE] = {0};
    int low_counts[OBSERVATION_STATE][OBSERVATION_STATE] = {0};
    int high_totals[OBSERVATION_STATE] = {0};
    int low_totals[OBSERVATION_STATE] = {0};

    for (int i = 0; i < seq_len - 1; i++) {
        int current_base = base_to_index(sequence[i]);
        int next_base = base_to_index(sequence[i + 1]);
        int current_state = hidden_states[i];

        if (current_base != -1 && next_base != -1) {
            if (current_state == 0) { // High 隱狀態
                high_counts[current_base][next_base]++;
                high_totals[current_base]++;
            } else if (current_state == 1) { // Low 隱狀態
                low_counts[current_base][next_base]++;
                low_totals[current_base]++;
            }
        }
    }

    for (int i = 0; i < OBSERVATION_STATE; i++) {
        if (high_totals[i] == 0) {
            printf("Warning: No transitions recorded for High state %d.\n", i);
        }
        for (int j = 0; j < OBSERVATION_STATE; j++) {
            high_matrix[i][j] = (high_totals[i] > 0) ? (double)high_counts[i][j] / high_totals[i] : 0;
        }
    }

    for (int i = 0; i < OBSERVATION_STATE; i++) {
        if (low_totals[i] == 0) {
            printf("Warning: No transitions recorded for Low state %d.\n", i);
        }
        for (int j = 0; j < OBSERVATION_STATE; j++) {
            low_matrix[i][j] = (low_totals[i] > 0) ? (double)low_counts[i][j] / low_totals[i] : 0;
        }
    }
}

void print_emission_matrix(double emission_matrix[HIDDEN_STATE][OBSERVATION_STATE]) {
    printf("\nEmission Matrix (Hidden * Observation):\n");
    printf("       A      C      G      T\n");
    for (int i = 0; i < HIDDEN_STATE; i++) {
        printf("%-5s", states[i]);
        for (int j = 0; j < OBSERVATION_STATE; j++) {
            printf("%6.2f ", emission_matrix[i][j]);
        }
        printf("\n");
    }
}
