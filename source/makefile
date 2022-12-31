all: main

main: main.o factorization.o 
	gcc -g -Wall -o main main.o factorization.o -lgmp

main.o : main.c
	gcc -g -Wall -c main.c -lgmp

factorization.o : factorization.c
	gcc -g -c -Wall factorization.c  -lgmp 

clean :
	rm -f  *.o *~