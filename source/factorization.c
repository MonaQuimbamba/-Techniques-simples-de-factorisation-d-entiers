
#include "Includes/factorization.h"


unsigned long int cputime()
{
    struct rusage rus;
    getrusage (0, &rus);
    return rus.ru_utime.tv_sec * 1000 + rus.ru_utime.tv_usec / 1000;
}

void print_primes_factors(PrimeFactors *f,mpz_t N){
    gmp_printf("%Zd = ", N);
     for (int i = 0; i < f->num_factors; i++){
                if (f->exponents[i] != 1) gmp_printf("%Zd^%d", f->prime_factors[i], f->exponents[i]);
                else gmp_printf("%Zd", f->prime_factors[i], f->exponents[i]);
                if (i != f->num_factors - 1) printf(" * ");
            }
            printf("\n\n");
   
}





bool p_minus_1(mpz_t n, mpz_t d, mpz_t B1, mpz_t B2){

 

    if (mpz_probab_prime_p(n, 10) > 0){ // first primality test 
         mpz_set(d, n);
        return true;
    }


    mpz_t a, p, q,tmp;
    gmp_randstate_t generateur;

    mpz_init(a);
    gmp_randinit_default(generateur);
    gmp_randseed_ui(generateur, time(NULL));
    mpz_urandomm(a, generateur, n);

    // Select a random integer a, 2 ≤ a ≤ n − 1, 
    if (mpz_cmp_ui(a, 1) <= 0) mpz_set_ui(a,2);

    //  compute d = gcd(a, n). If d ≥ 2
    //  then return(d).
    mpz_init(tmp);
    mpz_gcd(tmp, a, n);
    if (mpz_cmp_ui(tmp, 2) >= 0){
        mpz_set(d, tmp);
        mpz_clears(a, tmp, NULL);
        return true;
    }

    mpz_inits(p, q, NULL);
    mpz_set_ui(p, 2); // set the first prime number 
    mpz_set(d, n);
    bool found = false;
    bool go_next_stage=false;

    // stage 1 
    while((mpz_cmp(d, n) == 0) & (go_next_stage==false)) // make sure that if at the end d = n then we can do it again with a new base 'a'
    {
        while (mpz_cmp(p, B1) <= 0){ // make sure that p <= B1
            mpz_set_ui(q, 1);
            while(mpz_cmp(q, B1) <= 0){ // For each prime q ≤ B do the following:
                mpz_powm(a, a, p, n);
                mpz_mul(q, q, p); // Compute a←a^(q^l) mod n 
            }
            mpz_sub_ui(tmp, a, 1); // tmp <- a - 1
            mpz_gcd(d, tmp, n);    // d <- gcd(tmp,n)
            
            if (mpz_cmp_ui(d, 1) > 0){ // a divisor found then break 
                found = true;
                break;
            }else{ // otherwise  go to the next prime <= B1 and keep the same base 'a'      
              mpz_nextprime(p, p);
            }
           
        }
        if (mpz_cmp_ui(d, 1) == 0) go_next_stage = true;  // all prime powers are tested, no divisor found then go to stage 2 
        // d = n  then terminate the algorithm with failure try do it again with a new base 'a'
        else if (mpz_cmp(d, n) == 0){
            while (mpz_cmp_ui(a, 1) <= 0) mpz_urandomm(a, generateur, n);
            mpz_gcd(tmp, a, n);
            if (mpz_cmp_ui(tmp, 1) > 0){ // a divisor found then break 
                mpz_set(d, tmp);
                found = true;
                break;
            }
            mpz_set_ui(p, 2);
            mpz_set_ui(q, 1);
        } 
       
         
    }

    gmp_randclear(generateur);
    // stage  2 
    if ((mpz_cmp(B2, B1) > 0) && (go_next_stage == true)){
             
             mpz_t b,tmp_b,tmp_p,tmp_exp,t;
             mpz_inits(b,tmp_b,tmp_p,tmp_exp,t,NULL);
             mpz_set_ui(a,2);
             
            for(int i = 1 ; i <= mpz_get_ui(B1)  ; i++){ // computes ak+1 = ak^(k+1) mod n, for 1 <= k <= B1 − 1.
             mpz_powm_ui(a,a,i,n);
            }
                mpz_set(b,a); // set b(o)=a(B1)
                mpz_nextprime(p,B1);
                while(mpz_cmp(B2,p)>0){ // We then compute b1 = b(o)^l1  mod n, ...
                        mpz_powm(b,b,p,n);
                        mpz_sub_ui(tmp, b, 1);
                        mpz_gcd(tmp, tmp, n);
                        if (mpz_cmp_ui(tmp, 1) > 0 && mpz_cmp(tmp, n) < 0){
                            mpz_set(d, tmp);
                            found = true;
                            break;
                            
                        }
                       mpz_nextprime(p,p);
                }
                mpz_clear(b);

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
                
    }
    mpz_clears(a, p, q, tmp,NULL);
    return found;
}

int fact_p_1_pollard(mpz_t n,mpz_t B1,mpz_t B2,PrimeFactors *f)
{
   
  
    int res;
    mpz_t N, d, p, c;
    bool fact;
   
    mpz_init_set(N, n);

    mpz_inits(d, p, c, NULL);
    mpz_set(d, N);
    mpz_set_ui(p, 1);


    while (mpz_cmp_ui(N, 1) > 0){
        fact = p_minus_1(d, p, B1, B2);
        if (!fact) break;
       
        if (mpz_probab_prime_p(p, 10) > 0){  // a factor p is found then  test its primality 
               
               uint64_t exponent = 0;
                while (mpz_divisible_p(N, p) != 0) {
                mpz_divexact(N, N, p);
                exponent++;
                }                
                if (exponent > 0) {
                add_factor(f, p, exponent);
                }
             
            //add_prime_factor(f,N,p); // add the factor to the primes struct 
            mpz_set(d, N);           //  set d to the remaining factor N/p^i
            if (mpz_probab_prime_p(d, 10) > 0 || mpz_cmp_ui(d,1)==0){ //  primality test for d to check wheter we complete the factorization or not
                add_factor(f, d, exponent);
                mpz_set_ui(N,1);
                break;
            }
        }
        else mpz_set(d, p);  // in case factor p is not prime, apply Pollard's p-1 for p to find a prime factor of n
    }

 
    //  complete factorization
    if (mpz_cmp_ui(N, 1) == 0) res= 0;
    //  incomplete factorization or fail to factorized 
    else if (!fact)  res= -1 ; 
    
    
    mpz_clears(N, d, p, c, NULL);

    return res;
    
                
}


void add_factor(PrimeFactors* pf, mpz_t prime, uint64_t exponent) {
  pf->num_factors++;
  pf->prime_factors = (mpz_t*) realloc(pf->prime_factors, pf->num_factors * sizeof(mpz_t));
  pf->exponents = (uint64_t*) realloc(pf->exponents, pf->num_factors * sizeof(uint64_t));
  mpz_init_set(pf->prime_factors[pf->num_factors - 1], prime);
  pf->exponents[pf->num_factors - 1] = exponent;
}


bool trial_division(mpz_t p, mpz_t n, mpz_t p_max){
    if (mpz_probab_prime_p(n, 10) != 0){
        mpz_set(p,n);
        return true;
    }
    mpz_init_set_ui(p, 2); // set the first prime  we can next chose from which number we want to start 
    while (mpz_cmp(p, p_max) <= 0){
        if (mpz_divisible_p(n, p) != 0) {
            return true;
        }
        mpz_nextprime(p, p);
    }
    return false;
}

int fact_trialDivision(mpz_t n,mpz_t p_max,PrimeFactors *f)
{
   

   
    bool fact ;
    int res;
    mpz_t N, p;
    mpz_inits(p, N, NULL);
    mpz_set(N, n);

    while (mpz_cmp_ui(N, 1) > 0){
        // trial division by primes from p_min to p_max
        fact = trial_division(p, N,p_max);
        if (!fact) break; // end of factorizaction  or fail 
        if (mpz_probab_prime_p(p, 10) > 0){
               
                uint64_t exponent = 0;
                while (mpz_divisible_p(N, p) != 0) {
                mpz_divexact(N, N, p);
                exponent++;
                }
                if (exponent > 0) {
                add_factor(f, p, exponent);
                }
        }

    }
  
    //  complete factorization
    if (mpz_cmp_ui(N, 1) == 0) res= 0;
    //  incomplete factorization or fail to factorized 
    else if (!fact)  res= -1 ; 
    
    mpz_clears(N, p, NULL);

    return res;
    

}


bool pollard_rho_Floy_cycle(mpz_t n, mpz_t d,uint64_t nb_iterations){

    mpz_t t,x, y, c;
    mpz_inits(t,x, y, c, NULL);

    // Set the initial values for x and y
    mpz_set_ui(x, 2);
    mpz_set_ui(y, 2);

   
    mpz_set_ui(c,2); // x^2 + 2;
    

    // Set the initial value for d
    mpz_set_ui(d, 1);
    unsigned int i = 0;
    while(mpz_cmp_ui(d, 1) == 0){
        if (i > nb_iterations) break;
        i++;

        // Set x to (x^2 + c) mod n
        mpz_mul(x, x, x);
        mpz_add(x, x, c);
        mpz_mod(x, x, n);

        // Set y to (y^2 + c) mod n
        mpz_mul(y, y, y);
        mpz_add(y, y, c);
        mpz_mod(y, y, n);

        // Set y to (y^2 + c) mod n again
        mpz_mul(y, y, y);
        mpz_add(y, y, c);
        mpz_mod(y, y, n);
        
        // d = gcd(|x-y|,n)
        mpz_sub(d, x, y);
        mpz_abs(d, d);
        mpz_gcd(d, d, n);

        if ( mpz_cmp(d,n)!=0 && mpz_cmp_ui(d,1) !=0) {
            break;
        }
    }    
    mpz_clears(t,x, y, c, NULL);
    if ((mpz_cmp(d, n) == 0) || (mpz_cmp_ui(d, 1) == 0)) return false;
    else return true;
   
}      

bool pollard_rho_Brent_cycle(mpz_t n, mpz_t d, uint64_t nb_iterations){
    mpz_t x,y,c,t;
    mpz_inits(x,c,y,t,NULL);
   uint64_t r, k;

    // set random for y, c, m
    mpz_set_ui(c,1);

    // initialize y and p
    mpz_set_ui(y, 2);
    mpz_set_ui(d, 1);

    // r is a power of 2
    r = 1;
    
    while ((mpz_cmp_ui(d, 1) == 0) || (mpz_cmp(d, n) == 0)){
        
        // x to the current position of y, namely x_i
        mpz_set(x, y);

        k = 0;

        // y to x_(i+r)
        for (int i=0;i++;i<r){
            mpz_mul(y, y, y);
            mpz_add(y, y, c);
            mpz_mod(y, y, n);
        }

        // test gcd(|x-y|, n); where y ranges from x_(i+r) to x_(i+2r)
        while ((k < r) && (mpz_cmp_ui(d, 1) == 0)){
            mpz_mul(y, y, y);
            mpz_add(y, y, c);
            mpz_mod(y, y, n);
            mpz_sub(t, x, y);
            mpz_abs(t, t);
            mpz_mod(t, t, n);                  
            mpz_gcd(d, t, n);
            k++;
        }

        // going to the next power of 2
        r = r*2;
        if (r >= nb_iterations) break;
        if ( mpz_cmp(d,n)!=0 && mpz_cmp_ui(d,1) !=0) {
            break;
        }
    }

    
    mpz_clears(x,c,y,t,NULL);

    if ((mpz_cmp(d, n) == 0) || (mpz_cmp_ui(d, 1) == 0)) return false;
    else return true;
}




int fact_pollard_rho(mpz_t n,PrimeFactors *f,uint64_t nb_iterations){


    bool fact ;
    int res;
    mpz_t N, d, p;
    mpz_inits(d, p,N, NULL);
    mpz_set(N, n);
    mpz_set(d, n);
    mpz_set_ui(p, 1);

    while (mpz_cmp_ui(N, 1) > 0){
      
        fact = pollard_rho_Floy_cycle( d,p,  nb_iterations);
        // fail for both choices of polynomial
        if (!fact) break;

        // a factor p is found, perform primality test
        if (mpz_probab_prime_p(p, 10) != 0) {
            
                uint64_t exponent = 0;
                while (mpz_divisible_p(N, p) != 0) {
                    mpz_divexact(N, N, p);
                    exponent++;
                }                
                if (exponent > 0) add_factor(f, p, exponent);
                    
                mpz_set(d, N);

            // primality test for d
            if (mpz_probab_prime_p(d, 10) != 0){
                 add_factor(f, d, exponent);
                 mpz_set_ui(N,1);
                break;
            }
        } else mpz_set(d, p); // in case factor p is not prime, apply Pollard's rho for p to find a prime factor of n
    }


    // remaining factor is 1, complete factorization
    if (mpz_cmp_ui(N, 1) == 0) res = 0;
    // remaining factor is nontrivial and not prime, incomplete factorization
    else if ((!fact)) res = -1;

    mpz_clears(N, d, p, NULL);
    
    return res;


}