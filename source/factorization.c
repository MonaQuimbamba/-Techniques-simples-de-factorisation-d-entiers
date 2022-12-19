
#include "Includes/factorization.h"


unsigned long int cputime()
{
    struct rusage rus;
    getrusage (0, &rus);
    return rus.ru_utime.tv_sec * 1000 + rus.ru_utime.tv_usec / 1000;
}

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

void print_primes_factors(factor *f,mpz_t N){
    gmp_printf("%Zd = ", N);
     for (int i = 0; i < f->nb_factors; i++){
                if (f->exp[i] != 1) gmp_printf("%Zd^%d", f->prime_factors[i], f->exp[i]);
                else gmp_printf("%Zd", f->prime_factors[i], f->exp[i]);
                if (i != f->nb_factors - 1) printf(" * ");
            }
            printf("\n\n");
    // clear up all
     for (int i = 0; i < f->nb_factors; i++) mpz_clear(f->prime_factors[i]);
        free(f->prime_factors);
        free(f->exp);
}

void add_prime_factor(factor *f,mpz_t N,mpz_t p){

            // store the prime factor
           // uint64_t exp_N = f->exp[f->nb_factors-1];
            mpz_set(f->prime_factors[f->nb_factors-1], p);
            f->exp[f->nb_factors-1] = 0;

            // look how many times this primes number is prime factor of N 
            while(mpz_divisible_p(N, p)!=0){
                mpz_divexact(N, N, p);
                f->exp[f->nb_factors-1]++;
            }

           // gmp_printf(" =>  %Zu\n",N);
           // f->exp[f->nb_factors-1] *= exp_N;
            f->nb_factors++; // increase the nb  primes
            // realloc memory for the next factor
            f->prime_factors = (mpz_t*) realloc(f->prime_factors, (f->nb_factors)*(sizeof(mpz_t)));
            f->exp = (uint64_t*) realloc(f->exp, (f->nb_factors)*(sizeof(uint64_t)));
            // store the remaining factor N/p^e on the prime_factors in case this one is a prime number then we dont need to go back to algorithm to compute it
            mpz_init_set(f->prime_factors[f->nb_factors-1], N);
            f->exp[f->nb_factors-1] = 1;
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

void fact_p_1_pollard(mpz_t n,mpz_t B1,mpz_t B2,factor *f)
{


    
    printf("************* p-1 pollard ***************** \n");
    if (f->prime_factors == NULL)  f->prime_factors = (mpz_t*) malloc(sizeof(mpz_t));
    if (f->exp == NULL )  f->exp = (uint64_t*) malloc(sizeof(uint64_t));
    f->exp[0] = 1;
    f->nb_factors = 1;
    
    mpz_t N, d, p, c;
    mpz_init_set(N, n);

    mpz_inits(d, p, c, NULL);
    mpz_set(d, N);
    mpz_set_ui(p, 1);


    bool fact;
    while (mpz_cmp_ui(N, 1) > 0){
        fact = p_minus_1(d, p, B1, B2);
        if (!fact) break;
       
        if (mpz_probab_prime_p(p, 10) > 0){  // a factor p is found then  test its primality 
            add_prime_factor(f,N,p); // add the factor to the primes struct 
            mpz_set(d, N);           //  set d to the remaining factor N/p^i
            
            if (mpz_probab_prime_p(d, 10) > 0 || mpz_cmp_ui(d,1)==0){ //  primality test for d to check wheter we complete the factorization or not
                mpz_set_ui(N,1);
                break;
            }
        }
        else mpz_set(d, p);  // in case factor p is not prime, apply Pollard's p-1 for p to find a prime factor of n
    }

 
    // remaining factor is 1, complete factorization
    if (mpz_cmp_ui(N, 1) == 0){ printf("factorization complete \n");  print_primes_factors(f,n); }
    // remaining factor is nontrivial and not prime, incomplete factorization
    else if (!fact) {printf("factorization incomplete  try again with a upper bound B2 \n");  //print_primes_factors(f,N);
    }
    
    mpz_clears(N, d, p, c, NULL);
                
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

void fact_trialDivision(mpz_t n,mpz_t p_max,factor *f)
{
   
    if (f->prime_factors == NULL)  f->prime_factors = (mpz_t*) malloc(sizeof(mpz_t));
    if (f->exp == NULL )  f->exp = (uint64_t*) malloc(sizeof(uint64_t));
    f->exp[0] = 1;
    f->nb_factors = 1;
    bool fact ;
    mpz_t N, p;
    mpz_inits(p, N, NULL);
    mpz_set(N, n);

    while (mpz_cmp_ui(N, 1) > 0){
        // trial division by primes from p_min to p_max
        fact = trial_division(p, N,p_max);
        if (!fact) break; // end of factorizaction  or fail 
        add_prime_factor(f,N,p);
    }
    print_primes_factors(f,n);
    mpz_clears(N, p, NULL);
 
    

}
