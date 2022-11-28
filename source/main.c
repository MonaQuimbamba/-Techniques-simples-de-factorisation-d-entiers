
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gmp.h"
#include "Includes/p_1_pollard.h"
#include "Includes/rho_pollard.h"
#include "Includes/trial_division.h"


int main(int argc, char const *argv[])
{
    
    mpz_t N;
    mpz_init(N);
    mpz_set_str(N,"633564754957339397639948337059",0);
    fact_p_1_pollard(N);
    return 0;
}
