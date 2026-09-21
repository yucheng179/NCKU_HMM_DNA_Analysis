#ifndef MATRIX_UTILS_H
#define MATRIX_UTILS_H

#include "hmm.h"

void calculate_state_transition_matrix(char *sequence, int *hidden_states, int seq_len,
                                       double high_matrix[OBSERVATION_STATE][OBSERVATION_STATE],
                                       double low_matrix[OBSERVATION_STATE][OBSERVATION_STATE]);

void print_emission_matrix(double emission_matrix[HIDDEN_STATE][OBSERVATION_STATE]);

#endif
