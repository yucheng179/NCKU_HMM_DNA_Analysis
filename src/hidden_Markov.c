#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#define MAX_SEQ_LEN 1000000  // 假設最大序列長度
#define GEN_SEQ_LEN 100      // 生成的新序列長度
#define OBSERVATION_STATE 4  // 觀察符號數量 (A, C, G, T)
#define HIDDEN_STATE 2       // 隱狀態數量 (High CG, Low CG)

// 隱狀態名稱
const char *states[HIDDEN_STATE] = {"High", "Low"};
const char bases[OBSERVATION_STATE] = {'A', 'C', 'G', 'T'};

// 隱馬可夫模型參數
double initial_probabilities[HIDDEN_STATE] = {0.5, 0.5}; // 初始隱狀態機率
double transition_matrix[HIDDEN_STATE][HIDDEN_STATE] = {
    {0.8, 0.2},  // High CG -> {High CG, Low CG}
    {0.2, 0.8}   // Low CG -> {High CG, Low CG}
};
double emission_matrix[HIDDEN_STATE][OBSERVATION_STATE] = {
    {0.1, 0.4, 0.4, 0.1},  // High CG -> {A, C, G, T}
    {0.4, 0.1, 0.1, 0.4}   // Low CG -> {A, C, G, T}
};

// 將碱基字符映射到索引值
int base_to_index(char base) {
    switch (tolower(base)) {
        case 'a': return 0;
        case 'c': return 1;
        case 'g': return 2;
        case 't': return 3;
        default: return -1;  // 非法碱基
    }
}

// 隨機選擇一個狀態或符號
int random_choice(double *probabilities, int size) {
    double r = (double)rand() / RAND_MAX;
    double cumulative = 0;

    for (int i = 0; i < size; i++) {
        cumulative += probabilities[i];
        if (r < cumulative) {
            return i;
        }
    }
    return size - 1; // 預防累積誤差
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
            //printf("simplified lens++ len=%d\n", *len); // Debug: 打印長度
        }
    }

    fclose(file);

    if (*len == 0) {
        fprintf(stderr, "Error: The file %s contains no valid DNA sequence.\n", filename);
        free(sequence);
        exit(EXIT_FAILURE);
    }

    sequence[*len] = '\0'; // 添加字符串結束符
    return sequence;
}


void train_hmm(char *sequence, int seq_len, int *hidden_states) {
    int initial_counts[HIDDEN_STATE] = {0};
    int transition_counts[HIDDEN_STATE][HIDDEN_STATE] = {0};
    int emission_counts[HIDDEN_STATE][OBSERVATION_STATE] = {0};
    int state_counts[HIDDEN_STATE] = {0};

    // 計算初始機率和轉移/發射次數
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

    // 計算初始機率
    for (int i = 0; i < HIDDEN_STATE; i++) {
        initial_probabilities[i] = (double)initial_counts[i] / seq_len;
    }

    // 計算轉移機率
    for (int i = 0; i < HIDDEN_STATE; i++) {
        for (int j = 0; j < HIDDEN_STATE; j++) {
            transition_matrix[i][j] = (state_counts[i] > 0) ? 
                (double)transition_counts[i][j] / state_counts[i] : 0;
        }
    }

    // 計算發射機率
    for (int i = 0; i < HIDDEN_STATE; i++) {
        for (int j = 0; j < OBSERVATION_STATE; j++) {
            emission_matrix[i][j] = (state_counts[i] > 0) ? 
                (double)emission_counts[i][j] / state_counts[i] : 0;
        }
    }
}

// 生成隨機序列
void generate_sequence(int gen_len, char *observed_sequence, int *hidden_states) {
    int current_state = random_choice(initial_probabilities, HIDDEN_STATE);

    for (int i = 0; i < gen_len; i++) {
        // 保存隱狀態
        hidden_states[i] = current_state;

        // 根據隱狀態生成觀察符號
        int observed_base = random_choice(emission_matrix[current_state], OBSERVATION_STATE);
        observed_sequence[i] = bases[observed_base];

        // 根據隱狀態轉移到下一狀態
        current_state = random_choice(transition_matrix[current_state], HIDDEN_STATE);
    }

    observed_sequence[gen_len] = '\0'; // 添加字符串結束符
}

// 輸出生成序列中的隱狀態轉變
void print_generated_state_changes(int *hidden_states, int gen_len) {
    printf("\nGenerated State Changes:\n");
    for (int i = 1; i < gen_len; i++) {
        if (hidden_states[i] != hidden_states[i - 1]) {
            printf("State changed to %s at position %d\n", states[hidden_states[i]], i + 1);
        }
    }
}

// 生成特定序列的機率
double calculate_sequence_probability(char *sequence, int seq_len, int *hidden_states, 
                                      double high_matrix[OBSERVATION_STATE][OBSERVATION_STATE], 
                                      double low_matrix[OBSERVATION_STATE][OBSERVATION_STATE]) {

    double probability = initial_probabilities[hidden_states[0]] * 
                         emission_matrix[hidden_states[0]][base_to_index(sequence[0])];

    for (int t = 1; t < seq_len; t++) {
        int previous_state = hidden_states[t-1];
        int current_state = hidden_states[t];
        //int previous_base = base_to_index(sequence[t-1]);
        int current_base = base_to_index(sequence[t]);

        probability *= transition_matrix[previous_state][current_state];
        probability *= emission_matrix[current_state][current_base];
        
        /*if(current_state==0)
            probability *= high_matrix[previous_base][current_base];
        if(current_state==1)
            probability *= low_matrix[previous_base][current_base];*/


    }

    return probability;
}

// Viterbi 演算法
void viterbi(char *sequence, int seq_len, int *hidden_states) {
    double dp[seq_len][HIDDEN_STATE];
    int backtrack[seq_len][HIDDEN_STATE];

    // 初始化
    for (int i = 0; i < HIDDEN_STATE; i++) {
        dp[0][i] = initial_probabilities[i] * emission_matrix[i][base_to_index(sequence[0])];
        backtrack[0][i] = -1;
    }

    // 遞歸計算
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

    // 找出最終狀態
    double max_prob = -1;
    int last_state = -1;
    for (int i = 0; i < HIDDEN_STATE; i++) {
        if (dp[seq_len - 1][i] > max_prob) {
            max_prob = dp[seq_len - 1][i];
            last_state = i;
        }
    }

    // 回溯隱狀態序列
    hidden_states[seq_len - 1] = last_state;
    for (int t = seq_len - 2; t >= 0; t--) {
        hidden_states[t] = backtrack[t + 1][hidden_states[t + 1]];
    }
}

// 輸出隱狀態轉變信息
void print_hidden_state_changes(int *hidden_states, int seq_len) {
    printf("\nState Changes:\n");
    for (int i = 1; i < seq_len; i++) {
        if (hidden_states[i] != hidden_states[i - 1]) {
            printf("State changed to %s at position %d\n", states[hidden_states[i]], i + 1);
        }
    }
}

// 計算隱狀態內的轉移機率矩陣
void calculate_state_transition_matrix(char *sequence, int *hidden_states, int seq_len,
                                       double high_matrix[OBSERVATION_STATE][OBSERVATION_STATE],
                                       double low_matrix[OBSERVATION_STATE][OBSERVATION_STATE]) {
    int high_counts[OBSERVATION_STATE][OBSERVATION_STATE] = {0};
    int low_counts[OBSERVATION_STATE][OBSERVATION_STATE] = {0};
    int high_totals[OBSERVATION_STATE] = {0};
    int low_totals[OBSERVATION_STATE] = {0};

    for (int i = 0; i < seq_len - 1; i++) {
        int current_base = base_to_index(sequence[i]);
        int next_base = base_to_index(sequence[i + 1]);
        int current_state = hidden_states[i];

        if (current_base != -1 && next_base != -1) {
            if (current_state == 0) { // High 隱狀態
                high_counts[current_base][next_base]++;
                high_totals[current_base]++;
            } else if (current_state == 1) { // Low 隱狀態
                low_counts[current_base][next_base]++;
                low_totals[current_base]++;
            }
        }
    }

    // 計算 High 隱狀態的機率矩陣
    for (int i = 0; i < OBSERVATION_STATE; i++) {
        if (high_totals[i] == 0) {
            printf("Warning: No transitions recorded for High state %d.\n", i);
        }
        for (int j = 0; j < OBSERVATION_STATE; j++) {
            high_matrix[i][j] = (high_totals[i] > 0) ? (double)high_counts[i][j] / high_totals[i] : 0;
        }
    }

    // 計算 Low 隱狀態的機率矩陣
    for (int i = 0; i < OBSERVATION_STATE; i++) {
        if (low_totals[i] == 0) {
            printf("Warning: No transitions recorded for Low state %d.\n", i);
        }
        for (int j = 0; j < OBSERVATION_STATE; j++) {
            low_matrix[i][j] = (low_totals[i] > 0) ? (double)low_counts[i][j] / low_totals[i] : 0;
        }
    }

    

    
}
// 打印隱狀態發射矩陣
void print_emission_matrix(double emission_matrix[HIDDEN_STATE][OBSERVATION_STATE]) {
    printf("\nEmission Matrix (Hidden * Observation):\n");
    printf("       A      C      G      T\n"); // 碱基標籤
    for (int i = 0; i < HIDDEN_STATE; i++) {
        printf("%-5s", states[i]); // 隱狀態標籤
        for (int j = 0; j < OBSERVATION_STATE; j++) {
            printf("%6.2f ", emission_matrix[i][j]);
        }
        printf("\n");
    }
}


int main(int argc, char *argv[]) {
    srand(time(NULL)); // 初始化隨機數生成器

    // 訓練數據
    char *train_file = "./data/train_data.txt";
    //char *train_file = "./data/S288C_reference_sequence_R64-5-1_20240529.fsaa/S288C_reference_sequence_R64-5-1_20240529.fsa";
    int train_len;
    char *train_sequence = read_sequence_from_file(train_file, &train_len);

    // 測試數據
    char *test_file = "./data/vertified_data.txt";
    int test_len;
    char *test_sequence = read_sequence_from_file(test_file, &test_len);
    
    //printf("Simplified Sequence:\n%s\n", test_sequence);
    //printf("Sequence Length: %d\n", test_len);

    printf("\n\033[1;31mTraining Stage: \033[0m\n");
    printf("Observed Sequence:\n%s\n", train_sequence);

    int *hidden_states = malloc(train_len * sizeof(int));
    if (!hidden_states) {
        perror("Memory allocation failed");
        free(train_sequence);
        return EXIT_FAILURE;
    }

    // 使用 Viterbi 演算法計算隱狀態序列
    viterbi(train_sequence, train_len, hidden_states);

    // 輸出隱狀態轉變信息
    //print_hidden_state_changes(hidden_states, train_len);

    print_emission_matrix(emission_matrix);
    // 計算每個隱狀態的機率矩陣
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

     // 生成新序列
    /*char generated_sequence[GEN_SEQ_LEN + 1];
    int generated_states[GEN_SEQ_LEN];
    generate_sequence(GEN_SEQ_LEN, generated_sequence, generated_states);*/

    // 輸出生成的序列
    /*printf("\n\033[1;31mGenerating Stage: \033[0m\n");
    printf("Generated Sequence:\n%s\n", generated_sequence);*/

    // 輸出生成的隱狀態轉變
    //print_generated_state_changes(generated_states, GEN_SEQ_LEN);

     // 訓練階段
    int *train_states = malloc(train_len * sizeof(int));
    viterbi(train_sequence, train_len, train_states);
    train_hmm(train_sequence, train_len, train_states);
    free(train_states);

    //Debug
    printf("\nAfter Training: Transition Matrix:\n");
    for (int i = 0; i < HIDDEN_STATE; i++) {
        for (int j = 0; j < HIDDEN_STATE; j++) {
            printf("%.2f ", transition_matrix[i][j]);
        }
        printf("\n");
    }

    // 測試階段
    int *test_states = malloc(test_len * sizeof(int));
    viterbi(test_sequence, test_len, test_states);
    double test_probability = calculate_sequence_probability(test_sequence, test_len, test_states, high_matrix, low_matrix);

    // 輸出生成序列機率
    printf("\n");
    printf("\033[1;31mVertify Stage: \033[0m\n");
    printf("vertified Sequence:\n%s\n", test_sequence);
    printf("\n");

    print_hidden_state_changes(test_states, test_len);
    printf("\n");


    printf("Probability of model training: %.10e\n", test_probability);
    printf("Probability of random generation: %.10e\n", pow(4,test_len*(-1)));

    printf("\n");
    printf("improve rate: %g times\n", test_probability/(pow(4,test_len*(-1))));

    free(train_sequence);
    free(test_sequence);
    free(test_states);
    free(hidden_states);
    return 0;
}
