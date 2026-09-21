#include "hmm.h"
#include <stdio.h>
#include <ctype.h> 
#include <stdlib.h>
#include <math.h>

const char *states[HIDDEN_STATE] = {"High", "Low"};
const char bases[OBSERVATION_STATE] = {'A', 'C', 'G', 'T'};

double initial_probabilities[HIDDEN_STATE] = {0.5, 0.5};
double transition_matrix[HIDDEN_STATE][HIDDEN_STATE] = {
    {0.8, 0.2},
    {0.2, 0.8}
};
double emission_matrix[HIDDEN_STATE][OBSERVATION_STATE] = {
    {0.1, 0.4, 0.4, 0.1},
    {0.4, 0.1, 0.1, 0.4}
};

int base_to_index(char base) {
    switch (tolower(base)) {
        case 'a': return 0;
        case 'c': return 1;
        case 'g': return 2;
        case 't': return 3;
        default: return -1;
    }
}

int random_choice(double *probabilities, int size) {
    double r = (double)rand() / RAND_MAX;
    double cumulative = 0;

    for (int i = 0; i < size; i++) {
        cumulative += probabilities[i];
        if (r < cumulative) {
            return i;
        }
    }
    return size - 1;
}

void train_hmm(char *sequence, int seq_len, int *hidden_states) {
    int initial_counts[HIDDEN_STATE] = {0};
    int transition_counts[HIDDEN_STATE][HIDDEN_STATE] = {0};
    int emission_counts[HIDDEN_STATE][OBSERVATION_STATE] = {0};
    int state_counts[HIDDEN_STATE] = {0};

    initial_counts[hidden_states[0]]++;

    for (int i = 0; i < seq_len; i++) {
        int state = hidden_states[i];
        int observation = base_to_index(sequence[i]);

        if (observation != -1) {
            emission_counts[state][observation]++;
        }
        state_counts[state]++;

        if (i < seq_len - 1) {
            int next_state = hidden_states[i + 1];
            transition_counts[state][next_state]++;
        }
    }

    for (int i = 0; i < HIDDEN_STATE; i++) {
        initial_probabilities[i] = (double)initial_counts[i] / seq_len;
    }

    for (int i = 0; i < HIDDEN_STATE; i++) {
        for (int j = 0; j < HIDDEN_STATE; j++) {
            transition_matrix[i][j] = (state_counts[i] > 0) ? 
                (double)transition_counts[i][j] / state_counts[i] : 0;
        }
    }

    for (int i = 0; i < HIDDEN_STATE; i++) {
        for (int j = 0; j < OBSERVATION_STATE; j++) {
            emission_matrix[i][j] = (state_counts[i] > 0) ? 
                (double)emission_counts[i][j] / state_counts[i] : 0;
        }
    }
}

void viterbi(char *sequence, int seq_len, int *hidden_states) {
    double dp[seq_len][HIDDEN_STATE];
    int backtrack[seq_len][HIDDEN_STATE];

    for (int i = 0; i < HIDDEN_STATE; i++) {
        dp[0][i] = initial_probabilities[i] * emission_matrix[i][base_to_index(sequence[0])];
        backtrack[0][i] = -1;
    }

    for (int t = 1; t < seq_len; t++) {
        for (int s = 0; s < HIDDEN_STATE; s++) {
            double max_prob = -1;
            int prev_state = -1;

            for (int ps = 0; ps < HIDDEN_STATE; ps++) {
                double prob = dp[t - 1][ps] * transition_matrix[ps][s] * emission_matrix[s][base_to_index(sequence[t])];
                if (prob > max_prob) {
                    max_prob = prob;
                    prev_state = ps;
                }
            }

            dp[t][s] = max_prob;
            backtrack[t][s] = prev_state;
        }
    }

    double max_prob = -1;
    int last_state = -1;
    for (int i = 0; i < HIDDEN_STATE; i++) {
        if (dp[seq_len - 1][i] > max_prob) {
            max_prob = dp[seq_len - 1][i];
            last_state = i;
        }
    }

    hidden_states[seq_len - 1] = last_state;
    for (int t = seq_len - 2; t >= 0; t--) {
        hidden_states[t] = backtrack[t + 1][hidden_states[t + 1]];
    }
}

double calculate_sequence_probability(char *sequence, int seq_len, int *hidden_states, 
    double high_matrix[OBSERVATION_STATE][OBSERVATION_STATE], 
    double low_matrix[OBSERVATION_STATE][OBSERVATION_STATE]) {

    double probability = initial_probabilities[hidden_states[0]] * 
    emission_matrix[hidden_states[0]][base_to_index(sequence[0])];

    for (int t = 1; t < seq_len; t++) {
        int previous_state = hidden_states[t-1];
        int current_state = hidden_states[t];
        int current_base = base_to_index(sequence[t]);

        probability *= transition_matrix[previous_state][current_state];
        probability *= emission_matrix[current_state][current_base];
    }

    return probability;
}

void print_hidden_state_changes(int *hidden_states, int seq_len) {
    printf("\nState Changes:\n");
    for (int i = 1; i < seq_len; i++) {
        if (hidden_states[i] != hidden_states[i - 1]) {
            printf("State changed to %s at position %d\n", states[hidden_states[i]], i + 1);
        }
    }
}
