#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <inttypes.h>
#include <assert.h>
#include  <sys/types.h> 
#include  <sys/resource.h>
#include "gmp.h"


unsigned long int cputime();
void criblesimple(mpz_t * primes,uint64_t k);
void p_1_pollard(mpz_t N,mpz_t output);
/**
 * 
*/
void powerPrime(mpz_t pi,mpz_t ouput,uint64_t k);

/***
 *  return 2 if N is prime 
 *         1 if we found a the first prime factor  
 * */
int p_minus_1(mpz_t N,mpz_t output,mpz_t B1,mpz_t B2);
/**
 * 
*/
void fact_p_1_pollard(mpz_t N,uint16_t b1,uint64_t b2);
/**
 * 
*/
int fact_trialDivision(mpz_t *tab,mpz_t n,mpz_t pmax);

