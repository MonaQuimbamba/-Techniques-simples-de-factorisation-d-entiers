#include <stdio.h>
#include <gmp.h>

int main(int argc, char *argv[]) {
  // Initialize variables
  mpz_t n,p_max,p;
  mpz_inits(n,p,p_max,NULL);

  // Read in the number to be factorized from the command line
  mpz_set_str(n, argv[1], 10);
  mpz_set_str(p_max, argv[2], 10);                                                          ///   20   |  2
                                                                                           ///    10 
  mpz_init_set_ui(p, 2);

  

    //mpz_nextprime(p, p);
    while (mpz_cmp(p, p_max) <= 0){
        if (mpz_divisible_p(n, p) != 0) {
            gmp_printf(" p found = %Zu\n", p);
             mpz_divexact(n, n, p);
             //gmp_printf(" new n to test  = %Zu\n", n);
            
        }
        mpz_nextprime(p, p);
    }

  // Print the remaining factor
  gmp_printf("%Zd\n", n);

  // Clear memory and exit
  mpz_clear(n);
  mpz_clear(p);
  return 0;
}
