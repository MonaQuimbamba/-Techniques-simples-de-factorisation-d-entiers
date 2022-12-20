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



/// @brief struct use to store all primes factors of N
/// @param prime_factors : is the list of factors
/// @param  exp  : is the list of  exp of the factor 
/// @param  nb_factors: is the number of primes factros that N has

typedef struct factor {
    mpz_t *prime_factors;
    uint64_t *exp;
    int nb_factors;
} factor;


/// @brief function use to compute CPU time 
unsigned long int cputime();

/// @brief    function that allow to display all primes factors 
/// @param f  the struct which is made to store all primes factors
/// @param N  the number which we want factorize
void print_primes_factors(factor *f,mpz_t N);
/// @brief  function use to add a new prime factor to our struct and also it set the N value to the remaing of N/p^e
/// @param f the struct that is storing all primes factors
/// @param N the number which we are factorizing
/// @param p  the new prime factor found  
void add_prime_factor(factor *f,mpz_t N,mpz_t p);

/// @brief  function use to do the unitary test 
/// @param f  the struct with all primes factors 
/// @param N  the number which we're factorizing 
void testUnitaire(factor *f, mpz_t N);

/// @brief funtion that compute the P - 1 algorithm
/// @param n the number which we're factorizing 
/// @param p the prime factor found 
/// @param B1 the bound lower
/// @param B2  the bound upper
/// @return  true if we found a prime factor 
bool p_minus_1(mpz_t n,mpz_t p,mpz_t B1,mpz_t B2);

/// @brief funtion that compute the P - 1 algorithm
/// @param n the number which we're factorizing 
/// @param f  the struct with all primes factors 
/// @param B1 the bound lower
/// @param B2  the bound upper
int fact_p_1_pollard(mpz_t n,mpz_t  B1,mpz_t B2,factor *f);
/// @brief 
/// @param n 
/// @param p_max 
/// @param f 
/// @return 
int fact_trialDivision(mpz_t n,mpz_t p_max,factor *f);

/// @brief 
/// @param n 
/// @param p_max 
/// @param p 
/// @return 
bool trial_division(mpz_t p, mpz_t n,mpz_t p_max);