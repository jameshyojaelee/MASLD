CC = gcc
CFLAGS = -Wall -O3

all: prog1 prog2 prog3 prog4

prog1: prog1.c
	$(CC) $(CFLAGS) -o prog1 prog1.c

prog2: prog2.c
	$(CC) $(CFLAGS) -mavx2 -o prog2 prog2.c

prog3: prog3.c
	$(CC) $(CFLAGS) -mavx2 -o prog3 prog3.c

prog4: prog4.c
	$(CC) $(CFLAGS) -mavx2 -o prog4 prog4.c

clean:
	rm -f prog1 prog2 prog3 prog4