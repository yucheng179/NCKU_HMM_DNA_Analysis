CC = gcc
CFLAGS = -Wall -O2

TARGETS = zero_order first_order second_order

all: $(TARGETS)

zero_order: zero_order_Markov.c
	$(CC) $(CFLAGS) -o $@ $<

first_order: first_order_Markov.c
	$(CC) $(CFLAGS) -o $@ $<

second_order: second_order_Markov.c
	$(CC) $(CFLAGS) -o $@ $<

run_zero:
	./zero_order

run_first:
	./first_order

run_second:
	./second_order

clean:
	rm -f zero_order.exe first_order.exe second_order.exe
	
