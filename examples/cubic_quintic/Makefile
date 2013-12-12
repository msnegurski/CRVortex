CC=gcc
FLAGS=-Wall -O3 -ffast-math
INCLUDE=-I../../include
LINK=-L../../lib -lfft

all: M0 Mdef powers

M0: M0.c
	$(CC) $(FLAGS) $(INCLUDE) $(LINK) $^ -o $@

Mdef: Mdef.c
	$(CC) $(FLAGS) $(INCLUDE) $(LINK) $^ -o $@

powers: powers.c
	$(CC) $(FLAGS) $(INCLUDE) $(LINK) $^ -o $@

clean:
	rm -f *.txt
	rm -f M0 Mdef powers

