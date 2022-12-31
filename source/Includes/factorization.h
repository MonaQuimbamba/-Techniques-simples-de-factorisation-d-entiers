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


#ifndef STRUCT_FACTOR_H
#define STRUCT_FACTOR_H

/// @brief struct use to store all primes factors of N
/// @param prime_factors : is the list of factors
/// @param  exponents  : is the list of  exp of the factor 
/// @param  num_factors: is the number of primes factros that N has

typedef struct {
  mpz_t* prime_factors;
  uint64_t* exponents;
  size_t num_factors;
} PrimeFactors;
#endif

/// @brief function use to compute CPU time 
unsigned long int cputime();

/// @brief    function that allow to display all primes factors 
/// @param f  the struct which is made to store all primes factors
/// @param N  the number which we want factorize
void print_primes_factors(PrimeFactors *PrimeFactors,mpz_t N);



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
int fact_p_1_pollard(mpz_t n,mpz_t  B1,mpz_t B2,PrimeFactors *PrimeFactors);


/// @brief 
/// @param n 
/// @param p_max 
/// @param f 
/// @return 
int fact_trialDivision(mpz_t n,mpz_t p_max,PrimeFactors *PrimeFactors);

/// @brief 
/// @param n 
/// @param p_max 
/// @param p 
/// @return 
bool trial_division(mpz_t p, mpz_t n,mpz_t p_max);
/// @brief  function use to add a new prime factor to our struct and also it set the N value to the remaing of N/p^e
/// @param f the struct that is storing all primes factors
/// @param N the number which we are factorizing
/// @param p  the new prime factor found  

void add_factor(PrimeFactors* pf, mpz_t prime, uint64_t exponent);

bool pollard_rho_Floy_cycle(mpz_t n, mpz_t d,uint64_t nb_iterations);

bool pollard_rho_Brent_cycle(mpz_t n, mpz_t d, uint64_t nb_iterations);

int fact_pollard_rho(mpz_t n,PrimeFactors *f,uint64_t nb_iterations);