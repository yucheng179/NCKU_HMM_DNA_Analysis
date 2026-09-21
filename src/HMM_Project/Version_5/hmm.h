#ifndef HMM_H
#define HMM_H

#define OBSERVATION_STATE 4  
#define HIDDEN_STATE 4   


extern const char *states[HIDDEN_STATE];   
extern const char bases[OBSERVATION_STATE]; 

extern double initial_probabilities[HIDDEN_STATE];
extern double transition_matrix[HIDDEN_STATE][HIDDEN_STATE];
extern double emission_matrix[HIDDEN_STATE][OBSERVATION_STATE];
extern double initial_prob[HIDDEN_STATE];

int base_to_index(char base);
int random_choice(double *probabilities, int size);
double logsumexp(double a, double b);

void insert_forced_alu_segments(char *sequence, int seq_len, int *hidden_states);
void train_hmm(char *sequence, int seq_len, int *hidden_states);
void viterbi(char *sequence, int seq_len, int *hidden_states);

void baum_welch_train(char *sequence, int seq_len, int max_iter);
void forward_algorithm(char *sequence, int seq_len, double alpha[seq_len][HIDDEN_STATE]);
void forward_algorithm_with_mask(char *sequence, int seq_len, double alpha[seq_len][HIDDEN_STATE], int *hidden_states);
void backward_algorithm(char *sequence, int seq_len, double beta[seq_len][HIDDEN_STATE]);
void compute_gamma_xi(char *sequence, int seq_len, double alpha[seq_len][HIDDEN_STATE], double beta[seq_len][HIDDEN_STATE],
    double gamma[seq_len][HIDDEN_STATE], double xi[seq_len - 1][HIDDEN_STATE][HIDDEN_STATE]);
void baum_welch_train(char *sequence, int seq_len, int max_iter);

double calculate_sequence_probability(char *sequence, int seq_len, int *hidden_states,
                                      double high_matrix[OBSERVATION_STATE][OBSERVATION_STATE], 
                                      double low_matrix[OBSERVATION_STATE][OBSERVATION_STATE]);
void print_hidden_state_changes(int *hidden_states, int seq_len);

double compute_log_likelihood(char *sequence, int seq_len);
double compute_path_log_probability(char *sequence, int seq_len, int *states);
double forward_log_prob(char *sequence, int seq_len, int *hidden_states);
void print_state_usage_summary(int *states_seq, int seq_len);
void print_transition_entropy();
void print_exp_with_sign_and_zeros(int exponent);



#endif
