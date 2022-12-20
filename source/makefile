all: main

main: main.o factorization.o 
	gcc -Wall -o main main.o factorization.o -lgmp -lm

main.o : main.c
	gcc -Wall -c main.c -lgmp -lm

factorization.o : factorization.c
	gcc -c -Wall factorization.c -lm -lgmp

clean :
	rm -f  *.o *~