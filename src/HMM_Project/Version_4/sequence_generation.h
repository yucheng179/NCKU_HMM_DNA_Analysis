#ifndef SEQUENCE_GENERATION_H
#define SEQUENCE_GENERATION_H

#include "hmm.h"

void generate_sequence(int gen_len, char *observed_sequence, int *hidden_states);

#endif