#include "Includes/factorization.h"




int main(int argc, char const *argv[])
{

    
   if (argc < 4){
      printf("Use : %s n B1 B2  \n", argv[0]);
      exit(-1);
     }
    
    mpz_t N;
    mpz_t  B1;
    mpz_t B2;
    mpz_inits(N,B1,B2,NULL);
    
    mpz_set_str(N, argv[1], 10);
    mpz_set_str(B1, argv[2], 10);
    mpz_set_str(B2, argv[3], 10);

    factor f;
    f.prime_factors = NULL;
    f.exp = NULL;
    f.nb_factors = 0;

    fact_p_1_pollard(N,B1,B2,&f);
    //printf(" cputime %ld\n",cputime());
    testUnitaire(&f,N);
    print_primes_factors(&f,N);
    mpz_clears(N, B1, B2, NULL);

    
    

    return 0;
}
