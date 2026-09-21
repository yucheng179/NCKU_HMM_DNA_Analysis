#include "hmm.h"
#include <stdio.h>
#include <ctype.h> 
#include <stdlib.h>
#include <math.h>
#include <float.h>

#define EPSILON 1e-10

const char *states[HIDDEN_STATE] = {"High", "Low"};
const char bases[OBSERVATION_STATE] = {'A', 'C', 'G', 'T'};

double initial_probabilities[HIDDEN_STATE] = {0.5, 0.5};
double transition_matrix[HIDDEN_STATE][HIDDEN_STATE] = {
    {0.6, 0.4},
    {0.4, 0.6}
};
double emission_matrix[HIDDEN_STATE][OBSERVATION_STATE] = {
    {0.2, 0.3, 0.3, 0.2},
    {0.3, 0.2, 0.2, 0.3}
};


int base_to_index(char base) {
    switch (tolower(base)) {
        case 'a': return 0;
        case 'c': return 1;
        case 'g': return 2;
        case 't': return 3;
        default:
            fprintf(stderr, "[Error] Invalid base character: %c\n", base);
            exit(EXIT_FAILURE);
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

double logsumexp(double a, double b) {
    if (a == -INFINITY) return b;
    if (b == -INFINITY) return a;
    double max = (a > b) ? a : b;
    return max + log(exp(a - max) + exp(b - max));
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

    double alpha = 100.0;
    double min_transition = 0.01;

    for (int i = 0; i < HIDDEN_STATE; i++) {
        initial_probabilities[i] = (double)(initial_counts[i] + alpha) / (seq_len + alpha * HIDDEN_STATE);
    }

    for (int i = 0; i < HIDDEN_STATE; i++) {
        for (int j = 0; j < HIDDEN_STATE; j++) {
            transition_matrix[i][j] = (double)(transition_counts[i][j] + alpha) / 
                                      (state_counts[i] + alpha * HIDDEN_STATE);
            if (transition_matrix[i][j] < min_transition) {
                transition_matrix[i][j] = min_transition;
            }
        }
    }

    for (int i = 0; i < HIDDEN_STATE; i++) {
        for (int j = 0; j < OBSERVATION_STATE; j++) {
            emission_matrix[i][j] = (double)(emission_counts[i][j] + alpha) / 
                                    (state_counts[i] + alpha * OBSERVATION_STATE);
        }
    }
}


void viterbi(char *sequence, int seq_len, int *hidden_states) {
    double dp[seq_len][HIDDEN_STATE];
    int backtrack[seq_len][HIDDEN_STATE];

    for (int i = 0; i < HIDDEN_STATE; i++) {
        dp[0][i] = log(initial_probabilities[i]) + log(emission_matrix[i][base_to_index(sequence[0])]);
        backtrack[0][i] = -1;
    }

    for (int t = 1; t < seq_len; t++) {
        for (int s = 0; s < HIDDEN_STATE; s++) {
            double max_prob = -INFINITY;
            int prev_state = -1;

            for (int ps = 0; ps < HIDDEN_STATE; ps++) {
                double prob = dp[t - 1][ps] + log(transition_matrix[ps][s]) + log(emission_matrix[s][base_to_index(sequence[t])]);
                if (prob > max_prob) {
                    max_prob = prob;
                    prev_state = ps;
                }
            }

            dp[t][s] = max_prob;
            backtrack[t][s] = prev_state;
        }
    }

    double max_prob = -INFINITY;
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

void forward_algorithm(char *sequence, int seq_len, double alpha[seq_len][HIDDEN_STATE]) {
    for (int i = 0; i < HIDDEN_STATE; i++) {
        int obs_idx = base_to_index(sequence[0]);
        alpha[0][i] = log(initial_probabilities[i] + EPSILON) + log(emission_matrix[i][obs_idx] + EPSILON);
    }

    for (int t = 1; t < seq_len; t++) {
        int obs_idx = base_to_index(sequence[t]);
        for (int j = 0; j < HIDDEN_STATE; j++) {
            double sum_log = -INFINITY;
            for (int i = 0; i < HIDDEN_STATE; i++) {
                double log_prob = alpha[t - 1][i] + log(transition_matrix[i][j] + EPSILON);
                sum_log = logsumexp(sum_log, log_prob);
            }
            alpha[t][j] = sum_log + log(emission_matrix[j][obs_idx] + EPSILON);
        }
    }
}

void backward_algorithm(char *sequence, int seq_len, double beta[seq_len][HIDDEN_STATE]) {
    for (int i = 0; i < HIDDEN_STATE; i++) {
        beta[seq_len - 1][i] = 0;
    }

    for (int t = seq_len - 2; t >= 0; t--) {
        int obs_idx = base_to_index(sequence[t + 1]);
        for (int i = 0; i < HIDDEN_STATE; i++) {
            double sum_log = -INFINITY;
            for (int j = 0; j < HIDDEN_STATE; j++) {
                double log_prob = log(transition_matrix[i][j] + EPSILON) +
                                  log(emission_matrix[j][obs_idx] + EPSILON) +
                                  beta[t + 1][j];
                sum_log = logsumexp(sum_log, log_prob);
            }
            beta[t][i] = sum_log;
        }
    }
}

void compute_gamma_xi(char *sequence, int seq_len, double alpha[seq_len][HIDDEN_STATE], double beta[seq_len][HIDDEN_STATE],
                      double gamma[seq_len][HIDDEN_STATE], double xi[seq_len - 1][HIDDEN_STATE][HIDDEN_STATE]) {
    for (int t = 0; t < seq_len - 1; t++) {
        int next_obs = base_to_index(sequence[t + 1]);
        double denom = -INFINITY;

        for (int i = 0; i < HIDDEN_STATE; i++) {
            for (int j = 0; j < HIDDEN_STATE; j++) {
                double val = alpha[t][i] +
                             log(transition_matrix[i][j] + EPSILON) +
                             log(emission_matrix[j][next_obs] + EPSILON) +
                             beta[t + 1][j];
                denom = logsumexp(denom, val);
            }
        }

        for (int i = 0; i < HIDDEN_STATE; i++) {
            gamma[t][i] = -INFINITY;
            for (int j = 0; j < HIDDEN_STATE; j++) {
                xi[t][i][j] = alpha[t][i] +
                              log(transition_matrix[i][j] + EPSILON) +
                              log(emission_matrix[j][next_obs] + EPSILON) +
                              beta[t + 1][j] - denom;
                gamma[t][i] = logsumexp(gamma[t][i], xi[t][i][j]);
            }
        }
    }

    double denom_last = -INFINITY;
    for (int i = 0; i < HIDDEN_STATE; i++) {
        denom_last = logsumexp(denom_last, alpha[seq_len - 1][i]);
    }

    for (int i = 0; i < HIDDEN_STATE; i++) {
        gamma[seq_len - 1][i] = alpha[seq_len - 1][i] - denom_last;
    }
}

void baum_welch_train(char *sequence, int seq_len, int max_iter) {
    double (*alpha)[HIDDEN_STATE] = malloc(seq_len * sizeof(*alpha));
    double (*beta)[HIDDEN_STATE] = malloc(seq_len * sizeof(*beta));
    double (*gamma)[HIDDEN_STATE] = malloc(seq_len * sizeof(*gamma));
    double (*xi)[HIDDEN_STATE][HIDDEN_STATE] = malloc((seq_len - 1) * sizeof(*xi));

    if (!alpha || !beta || !gamma || !xi) {
        perror("Memory allocation failed for alpha/beta/gamma/xi");
        exit(EXIT_FAILURE);
    }

    double prev_log_likelihood = -INFINITY;

    for (int iter = 0; iter < max_iter; iter++) {
        forward_algorithm(sequence, seq_len, alpha);
        backward_algorithm(sequence, seq_len, beta);
        compute_gamma_xi(sequence, seq_len, alpha, beta, gamma, xi);

        // 計算 log-likelihood
        double log_likelihood = -INFINITY;
        for (int i = 0; i < HIDDEN_STATE; i++) {
            log_likelihood = logsumexp(log_likelihood, alpha[seq_len - 1][i]);
        }

        if (isnan(log_likelihood)) {
            fprintf(stderr, "[Error] NaN detected in log-likelihood at iteration %d\n", iter + 1);
            break;
        }

        printf("Iteration %2d: log-likelihood = %.6f (Δ=%.6f)\n",
               iter + 1, log_likelihood, log_likelihood - prev_log_likelihood);
        prev_log_likelihood = log_likelihood;

        // 更新 initial_probabilities
        for (int i = 0; i < HIDDEN_STATE; i++) {
            initial_probabilities[i] = exp(gamma[0][i]);
            if (isnan(initial_probabilities[i]) || initial_probabilities[i] < EPSILON) {
                //printf("⚠️  initial_probabilities[%d] = %.10f → set to EPSILON\n", i, initial_probabilities[i]);
                initial_probabilities[i] = EPSILON;
            }
        }

        // 更新 transition_matrix
        for (int i = 0; i < HIDDEN_STATE; i++) {
            double denom = -INFINITY;
            double min_transition = 0.01;
            for (int t = 0; t < seq_len - 1; t++)
                denom = logsumexp(denom, gamma[t][i]);

            for (int j = 0; j < HIDDEN_STATE; j++) {
                double numer = -INFINITY;
                for (int t = 0; t < seq_len - 1; t++)
                    numer = logsumexp(numer, xi[t][i][j]);

                transition_matrix[i][j] = exp(numer - denom);
                if (transition_matrix[i][j] < min_transition)
                    transition_matrix[i][j] = min_transition;
                if (isnan(transition_matrix[i][j]) || transition_matrix[i][j] < EPSILON) {
                    printf("⚠️  transition_matrix[%d][%d] invalid → set to EPSILON\n", i, j);
                    transition_matrix[i][j] = EPSILON;
                }
            }
            // ✅ Normalize row i
            double row_sum = 0.0;
            for (int j = 0; j < HIDDEN_STATE; j++)
                row_sum += transition_matrix[i][j];
            for (int j = 0; j < HIDDEN_STATE; j++)
                transition_matrix[i][j] /= row_sum;
        }

        // 更新 emission_matrix
        for (int i = 0; i < HIDDEN_STATE; i++) {
            double denom = -INFINITY;
            for (int t = 0; t < seq_len; t++)
                denom = logsumexp(denom, gamma[t][i]);

            for (int k = 0; k < OBSERVATION_STATE; k++) {
                double numer = -INFINITY;
                for (int t = 0; t < seq_len; t++) {
                    if (base_to_index(sequence[t]) == k)
                        numer = logsumexp(numer, gamma[t][i]);
                }

                emission_matrix[i][k] = exp(numer - denom);
                if (isnan(emission_matrix[i][k]) || emission_matrix[i][k] < EPSILON) {
                    printf("⚠️  emission_matrix[%d][%d] invalid → set to EPSILON\n", i, k);
                    emission_matrix[i][k] = EPSILON;
                }
            }
        }
    }

    free(alpha);
    free(beta);
    free(gamma);
    free(xi);
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
            //printf("State changed to %s at position %d\n", states[hidden_states[i]], i + 1);
        }
    }
}

double forward_log_prob(char *sequence, int seq_len) {
    double (*alpha)[HIDDEN_STATE] = malloc(seq_len * sizeof(*alpha));
    if (!alpha) {
        perror("Memory allocation failed for alpha");
        exit(EXIT_FAILURE);
    }

    forward_algorithm(sequence, seq_len, alpha);

    // log P(O) = logsumexp(alpha_T)
    double log_prob = -INFINITY;
    for (int i = 0; i < HIDDEN_STATE; i++) {
        log_prob = logsumexp(log_prob, alpha[seq_len - 1][i]);
    }

    free(alpha);
    return log_prob;
}

void print_state_usage_summary(int *states_seq, int seq_len) {
    int count[HIDDEN_STATE] = {0};

    for (int i = 0; i < seq_len; i++) {
        if (states_seq[i] >= 0 && states_seq[i] < HIDDEN_STATE) {
            count[states_seq[i]]++;
        }
    }

    printf("\n📊 Hidden State Usage Summary:\n");
    for (int i = 0; i < HIDDEN_STATE; i++) {
        double percentage = (100.0 * count[i]) / seq_len;
        printf("State %s (%d): %d times (%.2f%%)\n", states[i], i, count[i], percentage);
    }
}

void print_transition_entropy() {
    printf("\n🧮 Transition Matrix Entropy (per row):\n");

    for (int i = 0; i < HIDDEN_STATE; i++) {
        double entropy = 0.0;
        for (int j = 0; j < HIDDEN_STATE; j++) {
            double p = transition_matrix[i][j];
            if (p > 1e-10) {
                entropy -= p * log2(p);
            }
        }

        double max_entropy = log2(HIDDEN_STATE);
        printf("State %s (row %d): Entropy = %.4f / %.4f (%.1f%% randomness)\n",
               states[i], i, entropy, max_entropy, 100.0 * entropy / max_entropy);
    }
}

