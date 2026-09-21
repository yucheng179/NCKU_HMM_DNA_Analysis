#include "hmm.h"
#include <stdio.h>
#include <ctype.h> 
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

#define EPSILON 1e-10

#ifndef HAVE_STRNDUP
char *strndup(const char *s, size_t n) {
    char *result;
    size_t len = strnlen(s, n);  // 最多找 n 個長度
    result = malloc(len + 1);
    if (!result) return NULL;
    memcpy(result, s, len);
    result[len] = '\0';
    return result;
}
#endif

const char *states[HIDDEN_STATE] = {"High", "Low", "Mid", "Alu"};
const char bases[OBSERVATION_STATE] = {'A', 'C', 'G', 'T'};

double initial_probabilities[HIDDEN_STATE] = {0.33, 0.33, 0.34, 0.0};
double transition_matrix[HIDDEN_STATE][HIDDEN_STATE] = {
    { 0.6, 0.2, 0.2, 0.0 }, // High
    { 0.2, 0.6, 0.2, 0.0 }, // Low
    { 0.2, 0.2, 0.6, 0.0 }, // Mid
    { 0.3, 0.3, 0.4, 0.0 }  // Alu：可跳出，不能跳入
};
double emission_matrix[HIDDEN_STATE][OBSERVATION_STATE] = {
    {0.1, 0.4, 0.4, 0.1},   // High CG → 高 C/G
    {0.4, 0.1, 0.1, 0.4},   // Low CG → 高 A/T
    {0.25, 0.25, 0.25, 0.25}, // Mid → 均勻
    {0.27, 0.23, 0.23, 0.27}  // Alu → 可用固定序列訓練後改
};


int base_to_index(char base) {
    switch (tolower(base)) {
        case 'a': return 0;
        case 'c': return 1;
        case 'g': return 2;
        case 't': return 3;
        default:
            fprintf(stderr, "[Error] Invalid base character: %c (ASCII %d)\n", base, base);
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

int levenshtein_distance(const char *s1, const char *s2, int len, int max_dist) {
    int *v0 = malloc((len + 1) * sizeof(int));
    int *v1 = malloc((len + 1) * sizeof(int));
    for (int i = 0; i <= len; i++) v0[i] = i;

    for (int i = 0; i < len; i++) {
        v1[0] = i + 1;
        int min_row = v1[0];
        for (int j = 0; j < len; j++) {
            int cost = toupper(s1[i]) == toupper(s2[j]) ? 0 : 1;
            v1[j + 1] = fmin(fmin(v1[j] + 1, v0[j + 1] + 1), v0[j] + cost);
            if (v1[j + 1] < min_row) min_row = v1[j + 1];
        }
        if (min_row > max_dist) {
            free(v0); free(v1);
            return max_dist + 1;  // 立即放棄
        }
        int *tmp = v0; v0 = v1; v1 = tmp;
    }
    int result = v0[len];
    free(v0); free(v1);
    return result;
}


void insert_forced_alu_segments(char *sequence, int seq_len, int *hidden_states) {
    printf("[Alu Check] Start\n");

    const char *ALU_Y = "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGTGGATCACCTGAGGTCAGGAGTTCGAGACCAGCCTGGCCAACATGGTGAAACCCCGTCTCTACTAAAAATACAAAAATTAGCCGGGCGTGGTGGCGCGCGCCTGTAATCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATCGCTTGAACCCGGGAGGCGGAGGTTGCAGTGAGCCGAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGTGAGACTCCGTCTCAAAAAAAAAA";
    const char *ALU_S = "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGTGGATCACCTGAGGTCAGGAGTTCGAGACCAGCCTGGCCAACATGGTGAAACCCCGTCTCTACTAAAAATACAAAAATTAGCCGGGCGTGGTGGCGCGCGCCTGTAATCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATCGCTTGAACCCGGGAGGCGGAGGTTGCAGTGAGCCGAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAAAAAAA";
    const char *ALU_J = "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGTGGATCACCTGAGGTCAGGAGTTCGAGACCAGCCTGGCCAACATGGTGAAACCCCGTCTCTACTAAAAATACAAAAATTAGCCGGGCGTGGTGGCGCGCGCCTGTAATCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATCGCTTGAACCCGGGAGGCGGAGGTTGCAGTGAGCCGAGATCGTGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAAAAAAA";
    const char *ALUS[] = {ALU_Y, ALU_S, ALU_J};
    const int NUM_ALUS = 3;

    for (int k = 0; k < NUM_ALUS; k++) {
        const char *ALU = ALUS[k];
        int alu_len = strlen(ALU);
        int max_dist = (int)(alu_len * 0.07) + 1;

        for (int i = 0; i <= seq_len - alu_len; i++) {
            if (i % 10000 == 0) printf("🔍 Checking position %d/%d...\n", i, seq_len - alu_len);

            // 確保資料合法
            int invalid = 0;
            for (int j = 0; j < alu_len; j++) {
                char ch = toupper(sequence[i + j]);
                if (ch != 'A' && ch != 'C' && ch != 'G' && ch != 'T') {
                    invalid = 1;
                    break;
                }
            }
            if (invalid) continue;

            // 安全 substring
            char *window = strndup(&sequence[i], alu_len);
            int dist = levenshtein_distance(window, ALU, alu_len, max_dist);
            free(window);

            double mismatch_rate = (double)dist / (alu_len + 1e-6);
            if (mismatch_rate <= 0.07) {
                printf("⚠️  [%s] Alu-like segment at pos %d (%.2f%% mismatch)\n",
                       (k == 0 ? "AluY" : k == 1 ? "AluS" : "AluJ"), i, mismatch_rate * 100);
                for (int j = 0; j < alu_len; j++) {
                    hidden_states[i + j] = 3;
                }
                i += alu_len - 1;  // skip 避免重複偵測
            }
        }
    }

    printf("[Alu Check] End\n");
}


void train_hmm(char *sequence, int seq_len, int *hidden_states) {
    //printf("[Train HMM] Start\n");
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
    printf("[Train HMM] End\n");
}


void viterbi(char *sequence, int seq_len, int *hidden_states) {
    printf("[Viterbi] Start\n");

    // 🧠 使用 heap 分配避免 stack overflow
    double (*dp)[HIDDEN_STATE] = malloc(seq_len * sizeof(*dp));
    int (*backtrack)[HIDDEN_STATE] = malloc(seq_len * sizeof(*backtrack));

    if (!dp || !backtrack) {
        fprintf(stderr, "❌ [Error] Failed to allocate dp/backtrack\n");
        exit(EXIT_FAILURE);
    }

    int first_obs = base_to_index(sequence[0]);

    for (int i = 0; i < HIDDEN_STATE; i++) {
        dp[0][i] = log(initial_probabilities[i] + EPSILON) +
                   log(emission_matrix[i][first_obs] + EPSILON);
        backtrack[0][i] = -1;
    }

    for (int t = 1; t < seq_len; t++) {
        int obs = base_to_index(sequence[t]);

        for (int s = 0; s < HIDDEN_STATE; s++) {
            double max_prob = -INFINITY;
            int prev_state = -1;

            for (int ps = 0; ps < HIDDEN_STATE; ps++) {
                double trans = log(transition_matrix[ps][s] + EPSILON);
                double emit = log(emission_matrix[s][obs] + EPSILON);
                double prob = dp[t - 1][ps] + trans + emit;

                if (prob > max_prob) {
                    max_prob = prob;
                    prev_state = ps;
                }
            }

            if (prev_state == -1) {
                fprintf(stderr, "❌ [Error] No valid path to state %d at time %d\n", s, t);
                exit(EXIT_FAILURE);
            }

            dp[t][s] = max_prob;
            backtrack[t][s] = prev_state;

            if (isnan(dp[t][s]) || isinf(dp[t][s])) {
                fprintf(stderr, "❌ [Error] dp[%d][%d] is invalid: %.4f\n", t, s, dp[t][s]);
                exit(EXIT_FAILURE);
            }
        }
    }

    // 最後狀態找最大值
    double max_prob = -INFINITY;
    int last_state = -1;
    for (int i = 0; i < HIDDEN_STATE; i++) {
        if (dp[seq_len - 1][i] > max_prob) {
            max_prob = dp[seq_len - 1][i];
            last_state = i;
        }
    }

    if (last_state == -1) {
        fprintf(stderr, "❌ [Error] No valid final state at time %d\n", seq_len - 1);
        exit(EXIT_FAILURE);
    }

    // 回溯
    hidden_states[seq_len - 1] = last_state;
    for (int t = seq_len - 2; t >= 0; t--) {
        int next_state = hidden_states[t + 1];

        if (next_state < 0 || next_state >= HIDDEN_STATE) {
            fprintf(stderr, "❌ [Error] Invalid next_state: %d at time %d\n", next_state, t + 1);
            exit(EXIT_FAILURE);
        }

        hidden_states[t] = backtrack[t + 1][next_state];
        if (hidden_states[t] < 0 || hidden_states[t] >= HIDDEN_STATE) {
            fprintf(stderr, "❌ [Error] Backtrack returned invalid state: %d at time %d\n", hidden_states[t], t);
            exit(EXIT_FAILURE);
        }
    }

    free(dp);
    free(backtrack);

    printf("[Viterbi] End\n");
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

void forward_algorithm_with_mask(char *sequence, int seq_len, double alpha[seq_len][HIDDEN_STATE], int *hidden_states) {
    for (int i = 0; i < HIDDEN_STATE; i++) {
        int obs_idx = base_to_index(sequence[0]);
        alpha[0][i] = log(initial_probabilities[i] + EPSILON) + log(emission_matrix[i][obs_idx] + EPSILON);
    }

    for (int t = 1; t < seq_len; t++) {
        int obs_idx = base_to_index(sequence[t]);

        // 檢查是否為 Alu 區段
        if (hidden_states[t] == 3) {
            // 若為 Alu 區段，視為機率乘 1 → log 加 0（直接複製上一時刻的機率）
            for (int j = 0; j < HIDDEN_STATE; j++) {
                alpha[t][j] = alpha[t - 1][j];
            }
            continue;
        }
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
            if (i == 3) continue;  // 不更新 Alu 狀態轉出
            double denom = -INFINITY;
            double min_transition = 0.01;
            for (int t = 0; t < seq_len - 1; t++)
                denom = logsumexp(denom, gamma[t][i]);

            for (int j = 0; j < HIDDEN_STATE; j++) {
                if (j == 3) continue;  // 不更新任何轉入 Alu 的機率
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
            if (i == 3) continue;  // 不更新 Alu 的發射矩陣
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

// 新增：完整的 Forward-Backward 估計產生 P(O)
double compute_log_likelihood(char *sequence, int seq_len) {
    double (*alpha)[HIDDEN_STATE] = malloc(seq_len * sizeof(*alpha));
    forward_algorithm(sequence, seq_len, alpha);
    double log_likelihood = -INFINITY;
    for (int i = 0; i < HIDDEN_STATE; i++)
        log_likelihood = logsumexp(log_likelihood, alpha[seq_len - 1][i]);
    free(alpha);
    return log_likelihood;
}

// 新增：完整路徑機率的計算
// 透過解碼後的狀態序列 test_states 來還原 log P(O, Q | θ)
double compute_path_log_probability(char *sequence, int seq_len, int *states) {
    double log_prob = log(initial_probabilities[states[0]] + EPSILON) +
                      log(emission_matrix[states[0]][base_to_index(sequence[0])] + EPSILON);
    for (int t = 1; t < seq_len; t++) {
        int prev = states[t - 1];
        int curr = states[t];
        int obs = base_to_index(sequence[t]);
        log_prob += log(transition_matrix[prev][curr] + EPSILON) +
                    log(emission_matrix[curr][obs] + EPSILON);
    }
    return log_prob;
}

double forward_log_prob(char *sequence, int seq_len, int *hidden_states) {
    double (*alpha)[HIDDEN_STATE] = malloc(seq_len * sizeof(*alpha));
    if (!alpha) {
        perror("Memory allocation failed for alpha");
        exit(EXIT_FAILURE);
    }

    forward_algorithm_with_mask(sequence, seq_len, alpha, hidden_states);

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

// 自訂函式，輸出三位數正負號+補零格式指數字串
void print_exp_with_sign_and_zeros(int exponent) {
    char sign = (exponent >= 0) ? '+' : '-';
    int abs_exp = abs(exponent);
    printf("%c%03d", sign, abs_exp);
}

