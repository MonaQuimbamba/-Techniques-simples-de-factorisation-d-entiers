all: main

main: main.o trial_division.o rho_pollard.o p_1_pollard.o helper.o
	gcc -Wall -o main main.o trial_division.o rho_pollard.o  p_1_pollard.o helper.o

main.o : main.c
	gcc -Wall -c main.c -lgmp -lm

trial_division.o : trial_division.c
	gcc -c -Wall trial_division.c -lm -lgmp

rho_pollard.o : rho_pollard.c
	gcc -c -Wall rho_pollard.c -lm -lgmp

p_1_pollard.o : p_1_pollard.c
	gcc -c -Wall p_1_pollard.c -lm -lgmp 

helper.o : helper.c
	gcc -c -Wall helper.c -lm -lgmp 

clean :
	rm -f  *.o *~