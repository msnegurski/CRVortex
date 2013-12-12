CC=gcc
FLAGS=-Wall -O3 -ffast-math
INCLUDE=-I../../include
LINK=-L../../lib -lfft

all: dynamics

dynamics: dynamics.c
	$(CC)  $^ -o $@ $(FLAGS) $(INCLUDE) $(LINK)

clean:
	rm -f *.txt
	rm -f dynamics
