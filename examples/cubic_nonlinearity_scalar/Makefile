CC=gcc
FLAGS=-Wall -Ofast -fomit-frame-pointer -msse -msse2 -msse3 -mmmx -m3dnow
INCLUDE=-I../../include
LINK=-L../../lib -lfft

all: Mdef M0

Mdef: Mdef.c
	$(CC)  $^ -o $@ $(FLAGS) $(INCLUDE) $(LINK)

M0: M0.c
	$(CC)  $^ -o $@ $(FLAGS) $(INCLUDE) $(LINK)

clean:
	rm -f *.txt
	rm -f M0 Mdef

