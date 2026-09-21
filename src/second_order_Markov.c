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

// 計算 2 階起始狀態聯合機率 (兩碱基組合機率)
void calculate_initial_double_base_probabilities(char *sequence, int len, double init_probs[16]) {
    int counts[16] = {0};
    int total = 0;

    for (int i = 0; i < len - 1; i++) {
        int first = base_to_index(sequence[i]);
        int second = base_to_index(sequence[i + 1]);
        if (first != -1 && second != -1) {
            int state = first * 4 + second;
            counts[state]++;
            total++;
        }
    }

    for (int i = 0; i < 16; i++) {
        init_probs[i] = (total > 0) ? (double)counts[i] / total : 0;
    }
}

// 計算 2 階轉移矩陣
void calculate_second_order_matrix(char *sequence, int len, double matrix[16][4]) {
    int counts[16][4] = {0};  // 計算轉移次數
    int total[16] = {0};      // 計算以某個 2 碱基組合開頭的總數

    for (int i = 0; i < len - 2; i++) {
        int current1 = base_to_index(sequence[i]);
        int current2 = base_to_index(sequence[i + 1]);
        int next = base_to_index(sequence[i + 2]);

        if (current1 != -1 && current2 != -1 && next != -1) {
            int current_state = current1 * 4 + current2;
            counts[current_state][next]++;
            total[current_state]++;
        } else {
            printf("Warning: Invalid base detected at position %d, %d, %d\n", i, i + 1, i + 2);
        }
    }

    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < 4; j++) {
            matrix[i][j] = (total[i] > 0) ? (double)counts[i][j] / total[i] : 0;
        }
    }

    // 除錯訊息：顯示轉移矩陣
    char bases[4] = {'a', 'c', 'g', 't'};
    printf("Second-order Markov Transition Matrix:\n");
    for (int i = 0; i < 16; i++) {
        printf("State %c%c: ", bases[i / 4], bases[i % 4]);
        for (int j = 0; j < 4; j++) {
            printf("%c: %.2f ", bases[j], matrix[i][j]);
        }
        printf("\n");
    }
}

// 根據 2 階馬可夫模型隨機生成新序列
void generate_sequence_second_order(double matrix[16][4], int gen_len) {
    printf("Generated DNA sequence (2nd Order):\n");

    char bases[4] = {'a', 'c', 'g', 't'};
    int current1 = rand() % 4;  // 隨機選擇初始碱基
    int current2 = rand() % 4;

    printf("%c%c", bases[current1], bases[current2]);  // 打印前兩個碱基

    for (int i = 2; i < gen_len; i++) {
        int current_state = current1 * 4 + current2;
        double r = (double)rand() / RAND_MAX;
        double cumulative = 0;

        int next_base = -1;
        for (int j = 0; j < 4; j++) {
            cumulative += matrix[current_state][j];
            if (r < cumulative) {
                next_base = j;
                printf("%c", bases[j]);
                current1 = current2;
                current2 = j;
                break;
            }
        }

        if (next_base == -1) {
            printf("\nError: Failed to find next base at iteration %d. Random value: %.2f\n", i, r);
            break;
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

// 評估序列機率，加入起始雙碱基機率
void evaluate_probability_second(char *verify_seq, int verify_len, double matrix[16][4], double init_probs[16]) {
    if (verify_len < 2) return;

    double log_prob = 0;
    int first_state = base_to_index(verify_seq[0]) * 4 + base_to_index(verify_seq[1]);
    if (first_state >= 0 && first_state < 16 && init_probs[first_state] > 0) {
        log_prob += log2_custom(init_probs[first_state]);
    }

    int prev1 = base_to_index(verify_seq[0]);
    int prev2 = base_to_index(verify_seq[1]);
    for (int i = 2; i < verify_len; i++) {
        int curr = base_to_index(verify_seq[i]);
        if (prev1 != -1 && prev2 != -1 && curr != -1) {
            int state = prev1 * 4 + prev2;
            double prob = matrix[state][curr];
            if (prob > 0)
                log_prob += log2_custom(prob);
        }
        prev1 = prev2;
        prev2 = curr;
    }

    double log_rand = verify_len * log2_custom(0.25);
    print_probability_string("Probability (Second order Markov model)", log_prob);
    print_probability_string("Random generation probability (1/4^n)  ", log_rand);
    print_probability_string("Improvement rate                       ", log_prob - log_rand);
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

    // 計算 2 階轉移矩陣
    double second_order_matrix[16][4];
    calculate_second_order_matrix(train_sequence, train_len, second_order_matrix);

    // 計算起始雙碱基機率
    double init_probs[16];
    calculate_initial_double_base_probabilities(train_sequence, train_len, init_probs);

    // 初始化隨機數生成器
    srand(time(NULL));

    // 生成 2 階序列
    generate_sequence_second_order(second_order_matrix, GEN_SEQ_LEN);
    printf("\n");

    // 生成test序列的機率
    printf("Evaluating probability of the test sequence:\n");
    evaluate_probability_second(test_sequence, test_len, second_order_matrix, init_probs);
    printf("\n");

    free(train_sequence);
    free(test_sequence);
    return 0;
}
