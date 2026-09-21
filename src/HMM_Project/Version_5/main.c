#include "hmm.h"
#include "io.h"
#include "matrix_utils.h"
#include "sequence_generation.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main(int argc, char *argv[]) {
    srand(time(NULL)); // 初始化隨機數生成器

    // === Step 1: 讀取訓練與測試序列 ===
    char *train_file = "../../data/train_data.txt";
    int train_len;
    char *train_sequence = read_sequence_from_file(train_file, &train_len);

    char *test_file = "../../data/vertified_data.txt";
    int test_len;
    char *test_sequence = read_sequence_from_file(test_file, &test_len);

    // === Step 2: 初始模型估計（選擇性，用 Viterbi + train_hmm 快速初始化參數）===
    printf("\n\033[1;31m[Initial Parameter Estimation - Viterbi + train_hmm()]\033[0m\n");
    int *init_states = malloc(train_len * sizeof(int));
    viterbi(train_sequence, train_len, init_states);
    insert_forced_alu_segments(train_sequence, train_len, init_states);
    train_hmm(train_sequence, train_len, init_states);
    

    printf("\nInitial Transition Matrix:\n");
    for (int i = 0; i < HIDDEN_STATE; i++) {
        for (int j = 0; j < HIDDEN_STATE; j++) {
            printf("%.2f ", transition_matrix[i][j]);
        }
        printf("\n");
    }

    printf("\nInitial Emission Matrix:\n");
    print_emission_matrix(emission_matrix);

    // === Step 3: 使用 Baum-Welch EM 訓練參數 ===
    printf("\n\033[1;32m[Training Stage - Baum-Welch EM]\033[0m\n");
    baum_welch_train(train_sequence, train_len, 100);

    printf("\nTrained Transition Matrix (Baum-Welch):\n");
    for (int i = 0; i < HIDDEN_STATE; i++) {
        for (int j = 0; j < HIDDEN_STATE; j++) {
            printf("%.2f ", transition_matrix[i][j]);
        }
        printf("\n");
    }

    printf("\nTrained Emission Matrix (Baum-Welch):\n");
    print_emission_matrix(emission_matrix);

    // === Step 4: 對測試資料執行 Viterbi 解碼 ===
    int *test_states = malloc(test_len * sizeof(int));
    viterbi(test_sequence, test_len, test_states);

    // === Step 5: 計算測試資料的生成機率（log-prob）===
    double forward_logp = forward_log_prob(test_sequence, test_len, test_states);
    //double forward_total_prob = exp(forward_logp);
    /*double log_prob = log(initial_probabilities[test_states[0]]) + log(emission_matrix[test_states[0]][base_to_index(test_sequence[0])]);
    for (int i = 1; i < test_len; i++) {
        int prev = test_states[i - 1];
        int curr = test_states[i];
        int obs = base_to_index(test_sequence[i]);
        log_prob += log(transition_matrix[prev][curr]) + log(emission_matrix[curr][obs]);
    }
    double total_prob = exp(log_prob);*/

    // === Step 6: 比較與隨機產生的 baseline 機率 ===
    double log10_forward_prob = forward_logp / log(10);
    int exponent_forward = (int)floor(log10_forward_prob);
    double mantissa_forward = pow(10, log10_forward_prob - exponent_forward);

    double log_random_prob = -test_len * log(4);
    double log10_random_prob = log_random_prob / log(10);
    int exponent_random = (int)floor(log10_random_prob);
    double mantissa_random = pow(10, log10_random_prob - exponent_random);

    // 計算 improvement rate (用 log 機率相減，再轉成常數)
    double log_improvement = forward_logp - log_random_prob;
    double log10_improvement = log_improvement / log(10);
    int exponent_improvement = (int)floor(log10_improvement);
    double mantissa_improvement = pow(10, log10_improvement - exponent_improvement);

    printf("\n\033[1;35m[Verification and Probability Comparison]\033[0m\n");
    //printf("Probability (Baum-Welch trained model): %.10e\n", total_prob);
    //printf("Forward Probability (log):              %.6f\n", forward_logp);
    printf("FProbability (Baum-Welch trained model): %.10fe%+03d\n", mantissa_forward, exponent_forward);
    printf("Random generation probability (1/4^n):   %.10fe%+03d\n", mantissa_random, exponent_random);
    printf("Improvement rate:                        %.10fe", mantissa_improvement);
    print_exp_with_sign_and_zeros(exponent_improvement);
    printf("\n");

    /*printf("Forward Probability (log):             %.10e\n", forward_logp);
    printf("Random generation probability (log):  %.10e\n", log_random_prob);
    printf("Improvement rate (log):                %.10e\n", log_improvement);*/

    // === 額外分析：狀態使用比例與轉移熵 ===
    print_state_usage_summary(test_states, test_len);
    print_transition_entropy();

    // === 清理資源 ===
    free(init_states);
    free(train_sequence);
    free(test_sequence);
    free(test_states);

    return 0;
}
