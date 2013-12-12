CC=gcc
FLAGS=-Wall -Ofast -fomit-frame-pointer -msse -msse2 -msse3 -mmmx -m3dnow
INCLUDE=-I../include

all: libfft.a

libfft.a: fft_functions.o
	ar ruv $@ $^
	ranlib $@

fft_functions.o: fft_functions.c
	$(CC) $(FLAGS) $(INCLUDE) -c $^ -o $@

clean:
	rm *.o libfft.a
