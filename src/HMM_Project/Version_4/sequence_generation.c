#include "sequence_generation.h"
#include "hmm.h"
#include <stdio.h>
#include <stdlib.h>

void generate_sequence(int gen_len, char *observed_sequence, int *hidden_states) {
    int current_state = random_choice(initial_probabilities, HIDDEN_STATE);

    for (int i = 0; i < gen_len; i++) {
        hidden_states[i] = current_state;
        int observed_base = random_choice(emission_matrix[current_state], OBSERVATION_STATE);
        observed_sequence[i] = bases[observed_base];
        current_state = random_choice(transition_matrix[current_state], HIDDEN_STATE);
    }

    observed_sequence[gen_len] = '\0';
}