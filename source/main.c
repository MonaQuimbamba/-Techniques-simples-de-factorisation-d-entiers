#include "Includes/factorization.h"



int main(int argc, char const *argv[])
{
    
    mpz_t N;
    mpz_init(N);
    mpz_set_str(N,"1241143",0);
    fact_p_1_pollard(N);

	
    return 0;
}
