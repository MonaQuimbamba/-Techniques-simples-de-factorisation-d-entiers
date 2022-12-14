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

/**
 * 
 *  struct use to store all primes factors of N
 *  @prime_factors : is the list of factors
 *  @exp  : is the list of  exp of the factor 
 *  @ nb_factors: is the number of primes factros that N has
*/
typedef struct factor {
    mpz_t *prime_factors;
    unsigned long *exp;
    int nb_factors;
} factor;

/**
*   function use to compute CPU time 
*/
unsigned long int cputime();
/**
 * 
*/
void criblesimple(mpz_t * primes,uint64_t k);

/**
*/
void p_1_pollard(mpz_t N,mpz_t output);
/**
 * 
*/
void powerPrime(mpz_t pi,mpz_t ouput,uint64_t k);

/***
 *  return 2 if N is prime 
 *         1 if we found a the first prime factor  
 * */
int p_minus_1(mpz_t N,factor *fact,mpz_t B1,mpz_t B2);


/**
 * 
*/
int fact_p_m_1(mpz_t d, mpz_t n, mpz_t B1, mpz_t B2);

/**
 * 
*/
void fact_p_1_pollard(mpz_t N,uint16_t b1,uint64_t b2);
/**
 * 
*/
int fact_trialDivision(mpz_t *tab,mpz_t n,mpz_t pmax);

