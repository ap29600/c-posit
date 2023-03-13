CC    = gcc
FLAGS = -Wall -Wextra --pedantic -std=c11
OPT   = -O3 -march=native

test: main
	./main

main: main.c posit.h
	$(CC) -omain main.c $(FLAGS) $(OPT)

clean:
	$(RM) main
