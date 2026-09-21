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

    // 訓練數據
    char *train_file = "../../data/train_data.txt";
    int train_len;
    char *train_sequence = read_sequence_from_file(train_file, &train_len);

    // 測試數據
    char *test_file = "../../data/vertified_data.txt";
    int test_len;
    char *test_sequence = read_sequence_from_file(test_file, &test_len);


    printf("\n\033[1;31mTraining Stage: \033[0m\n");
    //printf("Observed Sequence:\n%s\n", train_sequence);

    // 分配記憶體存放隱狀態序列
    int *hidden_states = malloc(train_len * sizeof(int));
    if (!hidden_states) {
        perror("Memory allocation failed");
        free(train_sequence);
        return EXIT_FAILURE;
    }

    // 使用 Viterbi 演算法計算隱狀態序列
    viterbi(train_sequence, train_len, hidden_states);

    // 輸出 Transition矩陣
    printf("\nTransition Matrix:\n");
    for (int i = 0; i < HIDDEN_STATE; i++) {
        for (int j = 0; j < HIDDEN_STATE; j++) {
            printf("%.2f ", transition_matrix[i][j]);
        }
        printf("\n");
    }

    // 輸出隱狀態轉變信息
    print_emission_matrix(emission_matrix);

    // 計算每個隱狀態的轉移機率矩陣
    double high_matrix[OBSERVATION_STATE][OBSERVATION_STATE] = {0};
    double low_matrix[OBSERVATION_STATE][OBSERVATION_STATE] = {0};
    calculate_state_transition_matrix(train_sequence, hidden_states, train_len, high_matrix, low_matrix);

    // 輸出 High 隱狀態的機率矩陣
    printf("\nHigh State Transition Matrix:\n");
    for (int i = 0; i < OBSERVATION_STATE; i++) {
        for (int j = 0; j < OBSERVATION_STATE; j++) {
            printf("%.2f ", high_matrix[i][j]);
        }
        printf("\n");
    }

    // 輸出 Low 隱狀態的機率矩陣
    printf("\nLow State Transition Matrix:\n");
    for (int i = 0; i < OBSERVATION_STATE; i++) {
        for (int j = 0; j < OBSERVATION_STATE; j++) {
            printf("%.2f ", low_matrix[i][j]);
        }
        printf("\n");
    }

    // 訓練 HMM
    int *train_states = malloc(train_len * sizeof(int));
    if (!train_states) {
        perror("Memory allocation failed");
        free(train_sequence);
        free(hidden_states);
        return EXIT_FAILURE;
    }
    viterbi(train_sequence, train_len, train_states);
    train_hmm(train_sequence, train_len, train_states);
    free(train_states);

    int count_High = 0, count_Low = 0;
    for (int i = 0; i < train_len; i++) {
        if (hidden_states[i] == 0) count_High++;
        else count_Low++;
    }

    printf("\nDebug: Hidden State Distribution in Training Data:\n");
    printf("High: %d (%.2f%%), Low: %d (%.2f%%)\n", 
        count_High, (double)count_High / train_len * 100,
        count_Low, (double)count_Low / train_len * 100);


    // Debug: 輸出訓練後的轉移機率矩陣
    printf("\nAfter Training: Transition Matrix:\n");
    for (int i = 0; i < HIDDEN_STATE; i++) {
        for (int j = 0; j < HIDDEN_STATE; j++) {
            printf("%.2f ", transition_matrix[i][j]);
        }
        printf("\n");
    }

    // 測試階段
    int *test_states = malloc(test_len * sizeof(int));
    if (!test_states) {
        perror("Memory allocation failed");
        free(train_sequence);
        free(hidden_states);
        free(test_sequence);
        return EXIT_FAILURE;
    }
    viterbi(test_sequence, test_len, test_states);
    double test_probability = calculate_sequence_probability(test_sequence, test_len, test_states, high_matrix, low_matrix);

    

    // 輸出測試序列的機率
    printf("\n");
    printf("\033[1;31mVerification Stage: \033[0m\n");
    printf("Verified Sequence:\n%s\n", test_sequence);
    printf("\n");


    // Debug: 檢查 test_states 是否有變化
    printf("Debug: Hidden states for test sequence:\n");
    for (int i = 0; i < test_len; i++) {
        printf("%d ", test_states[i]);
    }
    printf("\n");
    
    print_hidden_state_changes(test_states, test_len);
    // Debug: 檢查 test_len 是否為 0
    printf("Debug: Test sequence length = %d\n", test_len);
    printf("\n");

    printf("Probability of model training: %.10e\n", test_probability);
    printf("Probability of random generation: %.10e\n", pow(4, test_len * (-1)));

    printf("\n");
    printf("Improvement rate: %g times\n", test_probability / pow(4, test_len * (-1)));

    // 釋放記憶體
    free(train_sequence);
    free(test_sequence);
    free(test_states);
    free(hidden_states);
    
    return 0;
}
