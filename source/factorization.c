
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
    mpz_t a, p, q,b, tmp;
    gmp_randstate_t generateur;

    mpz_init(a);
    gmp_randinit_default(generateur);
    gmp_randseed_ui(generateur, time(NULL));
    mpz_urandomm(a, generateur, n);

    // Select a random integer a, 2 ≤ a ≤ n − 1, 
    while (mpz_cmp_ui(a, 1) <= 0) mpz_urandomm(a, generateur, n);

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
             mpz_init(b);
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
                
    }
   
    // stage 2 of the algorithm, only started when B2 > B1 and no factor found after stage 1
   /* if ((mpz_cmp(B2, B1) > 0)){
        // precompute the prime gap table
        // gap between primes less than 10^15 is less than 1000

        mpz_t p_next, g, gap[1000];
        mpz_inits(p_next, g, NULL);

        // at the end of stage 1, a = base^(product of prime powers <= B1) mod n
        // compute a^(2k) for 0 < 2k < 1000 and store in gap[2k-1]
        for(int i = 1; i < 1000; i = i+2){
            mpz_init(gap[i]);
            mpz_powm_ui(gap[i], a, i+1, n);

        }
        // p is the first prime > B1, then first we must compute a^p mod n
        mpz_powm(a, a, p, n);


        // compute a = a^p * a^(next prime to p - p) for all primes p <= B2
        while(mpz_cmp(p, B2) <= 0){
            mpz_nextprime(p_next, p);
            mpz_sub(g, p_next, p); // compute gap = next prime p_next - current prime p
            int i = mpz_get_ui(g);
            mpz_mul(a, a, gap[i-1]); // compute a^p * a^gap
            mpz_mod(a, a, n);

            // compute gcd(a - 1, n)
            mpz_sub_ui(tmp, a, 1);
            mpz_gcd(d, tmp, n);
            if (mpz_cmp_ui(d, 1) > 0 && mpz_cmp(d, n) < 0){
                found = 0;
                break;
            }

            // go to the next prime
            mpz_set(p, p_next);
        }

        // clear data;
        mpz_clears(p_next, g, NULL);
        for(int i = 1; i < 1000; i = i+2) mpz_clear(gap[i]);
    }*/

    mpz_clears(a, p, q, tmp,NULL);
    return found;
}

void fact_p_1_pollard(mpz_t n,mpz_t B1,mpz_t B2,factor *f)
{

    // first primality test  in case N is prime
    if (mpz_probab_prime_p(n, 10) > 0){
            // N is already prime 
            gmp_printf(" is already prime ");

    } else{
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


    bool found ;
    while (mpz_cmp_ui(N, 1) > 0){
        found = p_minus_1(d, p, B1, B2);
        if (!found) break;
       
        if (mpz_probab_prime_p(p, 10) > 0){  // a factor p is found then  test its primality 
            add_prime_factor(f,N,p); // add the factor to the primes struct 
            mpz_set(d, N);           //  set d to the remaining factor N/p^i
            
            if (mpz_probab_prime_p(d, 10) > 0 || mpz_cmp_ui(d,1)==0){ //  primality test for d to check wheter we complete the factorization or not
                found = true;
                break;
            }
        }
        
        // else mpz_set(d, p);  // in case factor p is not prime, apply Pollard's p-1 for p to find a prime factor of n
    }

    // remaining factor is 1, complete factorization
    //if (mpz_cmp_ui(N, 1) == 0) r = 1;
    // remaining factor is nontrivial and not prime, incomplete factorization
    //else if ((r != 1) && (mpz_cmp(N, n) < 0)) r = 0;
    // else r = -1 and N = n indicate that algorithm fails to find a nontrivial factor of n -> failure

    mpz_clears(N, d, p, c, NULL);
                
    }     
    
   

    

}

int fact_trialDivision(mpz_t *tab,mpz_t n,mpz_t pmax)
{
   /* mpz_t r,p,quotient,result;
	mpz_inits(quotient,r,result,NULL);
	mpz_sqrt(r,n);
	mpz_init_set_ui(p,2);
	int e=0;
	int i=0;

	
	tab=(mpz_t*)(malloc(sizeof(mpz_t))); //allocation dynamique du tableau des facteurs de n
	int z=0;

	while(mpz_cmp(p,r)<=0)
	{
		mpz_mod(quotient,n,p);
		while(mpz_cmp_ui(quotient,0)==0)
		{
			mpz_cdiv_q(n,n,p); 
			e=e+1;
			mpz_mod(quotient,n,p);
		}
		if(e>0) // ie si p different de 1
		{
			mpz_init_set_ui(tab[i],1);
			mpz_pow_ui(result,p,e);
			if (mpz_cmp(result,pmax)<=0) // on teste si le facteur trouvé est inferieur ou pas à pmax
			{
				mpz_set(tab[i],result);
                
				mpz_sqrt(r,n);
				e=0;
				i+=1;
			}
		}
		mpz_add_ui(p,p,1); // on augmente p de 1
		
	}
	if(mpz_cmp_ui(n,1)>0 && mpz_cmp(n,pmax)<=0)
	{
		mpz_init_set_ui(tab[i],1);
		mpz_set(tab[i],n);
        
	}
	else
	{
        gmp_printf("%Zu ",tab[i]);
		i=i-1; // on ne retourne pas le dernier facteur de n qui est superieur à pmax
	}
	if (i!=0){ return 1;} else{	return -1;}*/
    return 0;

}
