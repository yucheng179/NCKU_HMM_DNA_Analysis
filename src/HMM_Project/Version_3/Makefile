CC = gcc
CFLAGS = -Wall -O2

all: hmm_project

hmm_project: main.o hmm.o io.o matrix_utils.o sequence_generation.o
	$(CC) $(CFLAGS) -o hmm_project main.o hmm.o io.o matrix_utils.o sequence_generation.o -lm

main.o: main.c hmm.h io.h matrix_utils.h sequence_generation.h
	$(CC) $(CFLAGS) -c main.c

hmm.o: hmm.c hmm.h
	$(CC) $(CFLAGS) -c hmm.c

io.o: io.c io.h
	$(CC) $(CFLAGS) -c io.c

matrix_utils.o: matrix_utils.c matrix_utils.h
	$(CC) $(CFLAGS) -c matrix_utils.c

sequence_generation.o: sequence_generation.c sequence_generation.h
	$(CC) $(CFLAGS) -c sequence_generation.c

clean:
	rm -f *.o hmm_project
