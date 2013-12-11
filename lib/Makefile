CC=gcc
FLAGS=-Wall -O3 -ffast-math
INCLUDE=-I../include

all: libfft.a

libfft.a: fft_functions.o
	ar ruv $@ $^
	ranlib $@

fft_functions.o: fft_functions.c
	$(CC) $(FLAGS) $(INCLUDE) -c $^ -o $@

clean:
	rm *.o libfft.a
