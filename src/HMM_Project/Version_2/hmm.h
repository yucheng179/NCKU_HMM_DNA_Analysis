#ifndef HMM_H
#define HMM_H

#define OBSERVATION_STATE 4  
#define HIDDEN_STATE 2       

extern const char *states[HIDDEN_STATE];   
extern const char bases[OBSERVATION_STATE]; 

extern double initial_probabilities[HIDDEN_STATE];
extern double transition_matrix[HIDDEN_STATE][HIDDEN_STATE];
extern double emission_matrix[HIDDEN_STATE][OBSERVATION_STATE];

int base_to_index(char base);
int random_choice(double *probabilities, int size);

void train_hmm(char *sequence, int seq_len, int *hidden_states);
void viterbi(char *sequence, int seq_len, int *hidden_states);
double calculate_sequence_probability(char *sequence, int seq_len, int *hidden_states,
                                      double high_matrix[OBSERVATION_STATE][OBSERVATION_STATE], 
                                      double low_matrix[OBSERVATION_STATE][OBSERVATION_STATE]);

void print_hidden_state_changes(int *hidden_states, int seq_len);

#endif
