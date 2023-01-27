#include "Includes/test_cases.h"

// gcc test.c -Wall -lgmp -o test 

int testUnitaire(factor *f, mpz_t n){  47446859972530264215606621291954070329666137695312
                                        115792089237316195423570985008687907853269984665640564039457584007913129639937
     int r=0;
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

//gmp_printf(" n= %Zu et n computed %Zu \n ",n,test);
if(mpz_cmp(test,n)!=0) r=-1;
mpz_clears(test,tmp,NULL);

return r;

}

int main(int argc, char const *argv[])
{
    mpz_t n,p_max,B1,B2;
    int res;
   int r;
   factor f;
   
    int mode = atoi(argv[1]);
    switch (mode)
   {
        case 1 :
            // usage  ./exec 1 n p_max   =>  for  Divisions successives
            mpz_inits(n,p_max,NULL);
            f.prime_factors = NULL;
            f.exp = NULL;
            f.nb_factors = 0;
            // input n, p
            mpz_set_str(n, argv[2], 10);
            mpz_set_str(p_max, argv[3], 10);
            r= fact_trialDivision(n,p_max,&f);
            if(r==0){
                res = testUnitaire(&f,n);
                if(res==0) printf("test passed ");
                else printf("test failled ");
            }  
            else if(r==-1)  printf("Factorization incompleted try again with a p_max  bigger \n");
            mpz_clears(n,p_max, NULL);
         
        break;

        case 2:
            // usage  ./exec 1 n p_max   =>  for  rho de Pollard
        break;

        case 3:
                    // usage  ./exec 1 n B1 B2   =>  for  p-1  de Pollard 
                    mpz_inits(n,B1,B2,NULL);
                    f.prime_factors = NULL;
                    f.exp = NULL;
                    f.nb_factors = 0;
                    // input n,B1,B2
                    mpz_set_str(n, argv[2], 10);
                    mpz_set_str(B1, argv[3], 10);
                    mpz_set_str(B2, argv[4], 10);
                    r = fact_p_1_pollard(n,B1,B2,&f);
                    if(r==0) {
                        res = testUnitaire(&f,n);
                        if(res==0) printf("test passed ");
                        else printf("test failled ");
                    } 
                    else if(r==-1)  printf("Factorization incompleted  try again with a bound B2 bigger \n");
                    mpz_clears(n, B1, B2, NULL);

                    
            break;

   }
    return res;
}

/* *********************  TO DO ****************************************
              mpz_set(b,a); // set b(o)=a(B1)
                mpz_nextprime(p,B1);
                int i=0;
                while(mpz_cmp(B2,p)>0){ // We then compute b1 = b(o)^l1  mod n, ...


                    if(i<1){ // to compute the b1 only 
                        

                        mpz_set(tmp_b,b);// set tmp_b <- b(k-1)
                        mpz_set(tmp_p,p); // set tmp_p <- p(k-1)

                        mpz_powm(b,b,p,n); // set b <- b(k+1)
                        mpz_sub_ui(tmp, b, 1);
                        mpz_gcd(tmp, tmp, n);
                        if (mpz_cmp_ui(tmp, 1) > 0 && mpz_cmp(tmp, n) < 0){
                            mpz_set(d, tmp);
                            found = true;
                            break;
                        }
                       mpz_nextprime(p,p);
                       i++;
                    }
                    else{
                          //bk+1 = bkclk+1−lk mod n
                          mpz_sub(tmp_exp,p,tmp_p);  // compute  expo  <- p(k+1) - p(k-1)
                          mpz_powm(t,tmp_b,tmp_exp,n); // t <- b(o)^(p - tmp_p) mod n 
                          mpz_mul(t,t,b); // t <- b(k+1) * b(k-1)
                          mpz_mod(t,t,n);  // b <- t mod n    => set b <- b(k+1)

                          mpz_set(tmp_b,b); // set tmp_b <- b(k - 1)
                          mpz_set(tmp_p,p); // set tmp_p <- p(k -1)
                          mpz_set(b,t);

                        mpz_sub_ui(tmp, b, 1);
                        mpz_gcd(tmp, tmp, n);
                        if (mpz_cmp_ui(tmp, 1) > 0 && mpz_cmp(tmp, n) < 0){
                            mpz_set(d, tmp);
                            found = true;
                            break;
                        }
                       mpz_nextprime(p,p);
                    }    
                }
                mpz_clears(b,tmp_b,tmp_p,tmp_exp,t,NULL);*/
                