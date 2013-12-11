CC=gcc
FLAGS=-Wall -O3 -ffast-math

all: front lambda_min lambda_max

front: front.c
	$(CC) $(FLAGS) $^ -o $@

lambda_min: lambda_min.c
	$(CC) $(FLAGS) $^ -o $@

lambda_max: lambda_max.c
	$(CC) $(FLAGS) $^ -o $@

clean:
	rm -f *.txt
	rm -f front lambda_min lambda_max
