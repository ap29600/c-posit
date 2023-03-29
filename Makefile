CC    = gcc
FLAGS = -Wall -Wextra --pedantic -std=c11 -lm -g3 # -fsanitize=undefined
OPT   = -O2

default: test

test: testbin
	./testbin

testbin: test.c testing_utils.h posit.h
	$(CC) -otestbin test.c $(FLAGS) $(OPT)

clean:
	$(RM) testbin
