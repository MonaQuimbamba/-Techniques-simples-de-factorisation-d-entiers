#include "Includes/factorization.h"



int main(int argc, char const *argv[])
{
    // 547 * 2269 = 1241143
    
    mpz_t N;
    mpz_init(N);
    uint64_t b1=7;
    uint64_t b2=13;
    mpz_set_str(N,"672038771836751227845696565342450315062141551559473564642434674541",0);
    fact_p_1_pollard(N,b1,b2);

   
  

    return 0;
}
