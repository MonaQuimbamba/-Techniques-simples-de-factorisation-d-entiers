#include "Includes/factorization.h"



void testUnitaire(factor *f, mpz_t N){
     mpz_t test,tmp;
     mpz_inits(test,tmp,NULL);
     mpz_set_ui(test,1);
     for (int i = 0; i < f->nb_factors; i++){
                 mpz_set_ui(tmp,1);
                if (f->exp[i] != 1) {  
                    mpz_pow_ui(tmp,f->prime_factors[i],f->exp[i]);
                }               
                else { 
                    mpz_set(tmp,f->prime_factors[i]);
                }  
                mpz_mul(test,test,tmp);
}

if(mpz_cmp(test,N)!=0) printf(" full factorization faill");
mpz_clears(test,tmp,NULL);

}