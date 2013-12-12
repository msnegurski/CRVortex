CC=gcc
FLAGS=-Wall -Ofast -fomit-frame-pointer -msse -msse2 -msse3 -mmmx -m3dnow
INCLUDE=-I../../include
LINK=-L../../lib -lfft

all: gauss_beam

gauss_beam: gauss_beam_collapse.c
	$(CC)  $^ -o $@ $(FLAGS) $(INCLUDE) $(LINK)

clean:
	rm -f *.txt
	rm -f gauss_beam
