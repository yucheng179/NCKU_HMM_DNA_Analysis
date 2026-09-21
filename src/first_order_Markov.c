#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#define MAX_SEQ_LEN 1000000  // 假設最大序列長度為 1,000,000
#define GEN_SEQ_LEN 100      // 生成的新序列長度

// 將碱基字符映射到索引值 (a -> 0, c -> 1, g -> 2, t -> 3)，非法字符返回 -1
int base_to_index(char base) {
    switch (tolower(base)) {
        case 'a': return 0;
        case 'c': return 1;
        case 'g': return 2;
        case 't': return 3;
        default: return -1;  // 非法碱基
    }
}

// 計算 1 階轉移矩陣
void calculate_first_order_matrix(char *sequence, int len, double matrix[4][4]) {
    int counts[4][4] = {0};  // 計算轉移次數
    int total[4] = {0};      // 計算以某個碱基開頭的總數

    for (int i = 0; i < len - 1; i++) {
        int current = base_to_index(sequence[i]);
        int next = base_to_index(sequence[i + 1]);
        if (current >= 0 && current < 4 && next >= 0 && next < 4) {
            counts[current][next]++;
            total[current]++;
        }
    }

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            matrix[i][j] = (total[i] > 0) ? (double)counts[i][j] / total[i] : 0;
        }
    }

    // 打印轉移矩陣
    char bases[4] = {'a', 'c', 'g', 't'};
    printf("First-order Markov Transition Matrix:\n");
    printf("    a      c      g      t\n");
    for (int i = 0; i < 4; i++) {
        printf("%c ", bases[i]);
        for (int j = 0; j < 4; j++) {
            printf("%6.2f ", matrix[i][j]);
        }
        printf("\n");
    }
}

// 計算起始碱基機率
void calculate_initial_base_probabilities(char *sequence, int len, double init_probs[4]) {
    int counts[4] = {0};
    int total = 0;

    for (int i = 0; i < len; i++) {
        int idx = base_to_index(sequence[i]);
        if (idx >= 0 && idx < 4) {
            counts[idx]++;
            total++;
        }
    }

    for (int i = 0; i < 4; i++) {
        init_probs[i] = (total > 0) ? (double)counts[i] / total : 0;
    }
}

// 根據 1 階馬可夫模型隨機生成新序列
void generate_sequence_first_order(double matrix[4][4], int gen_len) {
    printf("Generated DNA sequence (1st Order):\n");

    char bases[4] = {'a', 'c', 'g', 't'};
    int current = rand() % 4;  // 隨機選擇初始碱基

    for (int i = 0; i < gen_len; i++) {
        printf("%c", bases[current]);
        double r = (double)rand() / RAND_MAX;
        double cumulative = 0;

        for (int j = 0; j < 4; j++) {
            cumulative += matrix[current][j];
            if (r < cumulative) {
                current = j;
                break;
            }
        }
    }

    printf("\n");
}

void print_probability_string(const char *label, double log2_prob) {
    double log10_prob = log2_prob * log10(2);
    int exponent = (int)floor(log10_prob);
    double mantissa = pow(10, log10_prob - exponent);
    char sign = (exponent >= 0) ? '+' : '-';
    int abs_exp = abs(exponent);
    printf("%s: %.6fe%c%03d\n", label, mantissa, sign, abs_exp);
}

double log2_custom(double x) {
    return log(x) / log(2);
}

// 評估序列機率，加入起始碱基機率
void evaluate_probability_first(char *verify_seq, int verify_len, double matrix[4][4], double init_probs[4]) {
    if (verify_len == 0) return;

    double log_prob = 0;
    int first_base = base_to_index(verify_seq[0]);
    if (first_base >= 0 && init_probs[first_base] > 0) {
        log_prob += log2_custom(init_probs[first_base]);
    }

    int prev = first_base;
    for (int i = 1; i < verify_len; i++) {
        int curr = base_to_index(verify_seq[i]);
        if (prev != -1 && curr != -1) {
            double prob = matrix[prev][curr];
            if (prob > 0)
                log_prob += log2_custom(prob);
        }
        prev = curr;
    }

    double log_rand = verify_len * log2_custom(0.25);
    print_probability_string("Probability (First order Markov model)", log_prob);
    print_probability_string("Random generation probability (1/4^n) ", log_rand);
    print_probability_string("Improvement rate                      ", log_prob - log_rand);
}

// 從檔案讀取 DNA 序列
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
    return sequence;
}

int main(int argc, char *argv[]) {
    
    char *train_file = "./data/train_data.txt";
    int train_len;
    char *train_sequence = read_sequence_from_file(train_file, &train_len);

    char *test_file = "./data/vertified_data.txt";
    int test_len;
    char *test_sequence = read_sequence_from_file(test_file, &test_len);

    // 計算 1 階轉移矩陣
    double first_order_matrix[4][4];
    calculate_first_order_matrix(train_sequence, train_len, first_order_matrix);

    // 計算起始碱基機率
    double init_probs[4];
    calculate_initial_base_probabilities(train_sequence, train_len, init_probs);

    // 初始化隨機數生成器
    srand(time(NULL));

    // 生成 1 階序列
    generate_sequence_first_order(first_order_matrix, GEN_SEQ_LEN);
    printf("\n");

    // 生成test序列的機率
    printf("Evaluating probability of the test sequence:\n");
    evaluate_probability_first(test_sequence, test_len, first_order_matrix, init_probs);
    printf("\n");

    free(train_sequence);
    free(test_sequence);
    return 0;
}
