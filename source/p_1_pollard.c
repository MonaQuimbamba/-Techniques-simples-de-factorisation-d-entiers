
#include "Includes/p_1_pollard.h"



void p_1_pollard(mpz_t N,mpz_t output) 
{
     mpz_t d,a,minus_temp;
     mpz_inits(d,a,minus_temp,NULL);
     mpz_set_ui(a,2);
     mpz_set_ui(output,1);
     int iterations = 0;
     int j = 2;
    while (mpz_cmp_ui(N,j)>=0)
    {
       mpz_powm_ui(a, a,j,N);
       mpz_sub_ui(minus_temp,a,1);
       mpz_gcd (d, minus_temp, N);
       iterations++;
        if ( (mpz_cmp_ui(d,1) > 0) && (mpz_cmp(N,d) > 0)) {
            printf("nb iterations: %d\n", iterations);
            mpz_set(output,d);
            break;
        }
        j++;
    }
 
}

void fact_p_1_pollard(mpz_t N){
    mpz_t fact1,fact2;
    double duration;
    clock_t start, finish;
    mpz_inits(fact1,fact2,NULL);
    start = clock();
    p_1_pollard(N,fact1);
    finish = clock();
    mpz_divexact (fact2, N, fact1);
    gmp_printf("%Zu = %Zu * %Zu\n",N,fact1,fact2);

    duration = (double)(finish - start)/CLOCKS_PER_SEC;
    printf("Running time: %f seconds.\n\n", duration);
    
}
