#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <math.h>


#define MAX_SEQ_LEN 1000000  // 假設最大序列長度為 1,000,000
#define GEN_SEQ_LEN 100      // 生成的新序列長度

// 計算碱基出現概率的函數
void calculate_base_probabilities(char *sequence, int len, double *pA, double *pC, double *pG, double *pT) {
    int countA = 0, countC = 0, countG = 0, countT = 0;
    
    for (int i = 0; i < len; i++) {
        switch (tolower(sequence[i])) {
            case 'a': countA++; break;
            case 'c': countC++; break;
            case 'g': countG++; break;
            case 't': countT++; break;
        }
    }

    *pA = (double)countA / len;
    *pC = (double)countC / len;
    *pG = (double)countG / len;
    *pT = (double)countT / len;
}

// 根據碱基概率隨機生成新序列的函數
void generate_sequence(double pA, double pC, double pG, double pT, int gen_len) {
    printf("Generated DNA sequence (0th Order):\n");

    for (int i = 0; i < gen_len; i++) {
        double r = (double)rand() / RAND_MAX;

        if (r < pA) {
            printf("A");
        } else if (r < pA + pC) {
            printf("C");
        } else if (r < pA + pC + pG) {
            printf("G");
        } else {
            printf("T");
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

// 計算 log2(x)
double log2_custom(double x) {
    return log(x) / log(2);
}

// 評估零階馬可夫模型生成序列的機率
void evaluate_probability_zero(char *verify_seq, int verify_len, double pA, double pC, double pG, double pT) {
    double log_prob = 0;
    for (int i = 0; i < verify_len; i++) {
        switch (tolower(verify_seq[i])) {
            case 'a': log_prob += log2_custom(pA); break;
            case 'c': log_prob += log2_custom(pC); break;
            case 'g': log_prob += log2_custom(pG); break;
            case 't': log_prob += log2_custom(pT); break;
        }
    }
    double log_rand = verify_len * log2_custom(0.25);
    print_probability_string("Probability (Zero order Markov model)", log_prob);
    print_probability_string("Random generation probability (1/4^n)", log_rand);
    print_probability_string("Improvement rate                     ", log_prob - log_rand);
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


    // 計算各碱基的出現概率
    double pA, pC, pG, pT;
    calculate_base_probabilities(train_sequence, train_len, &pA, &pC, &pG, &pT);

    printf("Base probabilities:\n");
    printf("P(a) = %.4f\n", pA);
    printf("P(c) = %.4f\n", pC);
    printf("P(g) = %.4f\n", pG);
    printf("P(t) = %.4f\n", pT);

    // 初始化隨機數生成器
    srand(time(NULL));

    // 生成新的 DNA 序列
    generate_sequence(pA, pC, pG, pT, GEN_SEQ_LEN);
    printf("\n");

    // 生成test序列的機率
    printf("Evaluating probability of the test sequence:\n");
    //printf("Test sequence: %s\n", test_sequence);
    evaluate_probability_zero(test_sequence, test_len, pA, pC, pG, pT);
    printf("\n");

    free(train_sequence);
    free(test_sequence);
    return 0;
}
