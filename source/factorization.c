
#include "Includes/factorization.h"


unsigned long int cputime()
{
struct rusage rus;
getrusage (0, &rus);
return rus.ru_utime.tv_sec * 1000 + rus.ru_utime.tv_usec / 1000;
}

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


int fact_p_m_1(mpz_t d, mpz_t n, mpz_t B1, mpz_t B2)
{
    mpz_t a, p, q, t;
    mpz_init_set_ui(a, 1);

    gmp_randstate_t state;
    unsigned long int seed = time(NULL);
    gmp_randinit_default(state);
    gmp_randseed_ui(state, seed);
    mpz_urandomm(a, state, n);

    // randomize the base a
    while (mpz_cmp_ui(a, 1) <= 0) mpz_urandomm(a, state, n);
    
    // return gcd(a, n) if nontrivial
    mpz_init(t);
    mpz_gcd(t, a, n);
    if (mpz_cmp_ui(t, 1) > 0){
        mpz_set(d, t);
        mpz_clears(a, t, NULL);
        return 0;
    }

    mpz_inits(p, q, NULL);
    mpz_set_ui(p, 2);
    mpz_set(d, n);

    int r = -1;
    // stage 1
    while(mpz_cmp(d, n) == 0){
        while (mpz_cmp(p, B1) <= 0){

            mpz_set_ui(q, 1);
            while(mpz_cmp(q, B1) <= 0){
                mpz_powm(a, a, p, n); // a <- a^p mod n
                mpz_mul(q, q, p); // q <- q*p
            }

            mpz_sub_ui(t, a, 1);
            mpz_gcd(d, t, n);

            // a divisor > 1 found
            if (mpz_cmp_ui(d, 1) > 0){
                r = 0;
                break;
            }

            // go to the next prime
            mpz_nextprime(p, p);
        }

        // all prime powers are tested, no divisor found
        if (mpz_cmp_ui(d, 1) == 0) r = -1;

        // found divisor d = n, rerandomize the base a and restart the calculation
        else if (mpz_cmp(d, n) == 0){
            while (mpz_cmp_ui(a, 1) <= 0) mpz_urandomm(a, state, n);
            mpz_gcd(t, a, n);

            if (mpz_cmp_ui(t, 1) != 0){
                mpz_set(d, t);
                r = 0;
                break;
            }
            mpz_set_ui(p, 2);
            mpz_set_ui(q, 1);
        }
    }

    gmp_randclear(state);

    // stage 2 of the algorithm, only started when B2 > B1 and no factor found after stage 1
    if ((mpz_cmp(B2, B1) > 0) && (r != 0)){
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
            mpz_sub_ui(t, a, 1);
            mpz_gcd(d, t, n);
            if (mpz_cmp_ui(d, 1) > 0 && mpz_cmp(d, n) < 0){
                r = 0;
                break;
            }

            // go to the next prime
            mpz_set(p, p_next);
        }

        // clear data;
        mpz_clears(p_next, g, NULL);
        for(int i = 1; i < 1000; i = i+2) mpz_clear(gap[i]);
    }
    
    mpz_clears(a, p, q, t, NULL);
    return r;
}

int p_minus_1(mpz_t N,factor *fact,mpz_t B1,mpz_t B2)
{
   
    
    if (fact == NULL) fact = (factor*) malloc(sizeof(factor));
    fact->prime_factors = (mpz_t*) malloc(sizeof(mpz_t));
    fact->exp = (unsigned long*) malloc(sizeof(unsigned long));
    fact->exp[0] = 1;
    fact->nb_factors = 1;

    mpz_t d, p, c;
    mpz_inits(d, p, c, NULL);
    mpz_set(d, N);
    mpz_set_ui(p, 1);

    int r = 0;

    while (mpz_cmp_ui(N, 1) > 0){
        r = fact_p_m_1(p, d, B1, B2);
        if (r != 0) break;
        
        // a factor p is found, perform primality test
        if (mpz_probab_prime_p(p, 10) != 0){
            // store the prime factor
            unsigned long exp_N = fact->exp[fact->nb_factors-1];
            mpz_set(fact->prime_factors[fact->nb_factors-1], p);
            fact->exp[fact->nb_factors-1] = 0;
            // compute the exponent
            while(mpz_divisible_p(N, p)){
                mpz_divexact(N, N, p);
                fact->exp[fact->nb_factors-1]++;
            }

            fact->exp[fact->nb_factors-1] *= exp_N;

            // store the remaining factor N/p^e
            fact->nb_factors++;
            fact->prime_factors = (mpz_t*) realloc(fact->prime_factors, (fact->nb_factors)*(sizeof(mpz_t)));
            fact->exp = (unsigned long*) realloc(fact->exp, (fact->nb_factors)*(sizeof(unsigned long)));
            mpz_init_set(fact->prime_factors[fact->nb_factors-1], N);
            fact->exp[fact->nb_factors-1] = exp_N;

            // set d to the remaining factor N/p^e
            mpz_set(d, N);

            // primality test for d
            if (mpz_probab_prime_p(d, 10) != 0){
                // remaining factor is prime, hence complete factorization
                r = 1;
                break;  
            }          
        } else mpz_set(d, p);  // in case factor p is not prime, apply Pollard's p-1 for p to find a prime factor of n
    }

    // remaining factor is 1, complete factorization
    //if (mpz_cmp_ui(N, 1) == 0) r = 1;
    // remaining factor is nontrivial and not prime, incomplete factorization
    //else if ((r != 1) && (mpz_cmp(N, n) < 0)) r = 0;
    // else r = -1 and N = n indicate that algorithm fails to find a nontrivial factor of n -> failure 

    mpz_clears(N, d, p, c, NULL);

   /* mpz_t q,a,base,fact,tmp;
    mpz_t *primes;
    mpz_inits(q,a,base,fact,tmp,NULL);
    mpz_set_ui(q,1);
    //mpz_set_ui(a,2);

    gmp_randstate_t state;
    unsigned long int seed = time(NULL);
    gmp_randinit_default(state);
    gmp_randseed_ui(state, seed);
    mpz_urandomm(base, state, N);

    // randomize the base a
    // Select a random integer a, 2 ≤ a ≤ n − 1, 
    // to do !!!!!!!!
    while (mpz_cmp_ui(base, 1) <= 0) mpz_urandomm(base, state, N);
 
    //and compute d = gcd(a, n). If d ≥ 2 then return(d).
    mpz_gcd(tmp,base,N);
    if (mpz_cmp_ui(tmp,2)>=0) {
         mpz_set(output,tmp);
        return 2;   
    }

    
    // stage 1
    while(true)
    
    {   uint64_t b = mpz_get_ui(B1);
        primes = malloc(sizeof(mpz_t) * b-1);
        criblesimple(primes,b); // compute all primes lower than b 
        //for(int i =1 ; i < b; i++) gmp_printf(" %Zu ",primes[i]);
        int i=1; // compute the Q 
        // Let Q be the least common multiple of all powers of primes ≤ B that are ≤ n. 
        while(i < b) {
            if(mpz_cmp_ui(primes[i],0)!=0){          
                powerPrime(primes[i],q,b);
            }
            i++;
        }
        mpz_powm(a,base,q,N); // a <- a^q mod N
        mpz_sub_ui(a,a,1); // a <- a - 1
        mpz_mod(a,a,N); // a <- a mod N 
        mpz_powm_ui(fact,a,1,N); // fact <- a ^ 1 mod N 
        if(mpz_cmp_ui(fact,0)!=0){
            mpz_gcd(fact,a,N); // gcd(a,N)
            if(mpz_cmp_ui(fact,1)==0){
                // go to next  primes b 
                //mpz_set_ui(tmp,b); // tmp <- b
                mpz_nextprime(B1,B1); // tmp <- next prime
                //b= mpz_get_ui(tmp); // b <- next prime
                mpz_set_ui(q,1);
                
            }else{ //
                    
                    if(mpz_probab_prime_p(fact,10)>0){
                        mpz_set(output,fact);
                        return 1;
                    }
            }
            
        }
    }*/

    // stage 2 

}

void powerPrime(mpz_t pi,mpz_t ouput,uint64_t b )
{
    mpz_t tmp,test;
    mpz_inits(tmp,test,NULL);
    mpz_set(tmp,pi);
    int go=1;
    while(1){
        mpz_mul(test,tmp,pi);
        if(mpz_cmp_ui(test,b)!=-1){ mpz_set(tmp,tmp);    break;}
        else{mpz_set(tmp,test); go++;}
    }
    if(go==1) mpz_mul(ouput,ouput,pi);
     else mpz_mul(ouput,ouput,tmp);
}

void criblesimple(mpz_t * primes,uint64_t k){

    mpz_t testBoucle,mult,j,minus;
    mpz_inits(testBoucle,mult,j,minus,NULL);
     int  i;  
    for(i = 0; i < k; i++)  mpz_init_set_ui(primes[i],i+1);

    mpz_mul(mult,primes[1],primes[1]);
    mpz_set(testBoucle,mult);
    i=1;
    while( mpz_cmp_ui(testBoucle,k)<0)  
    {  
        if(mpz_cmp_ui(primes[i],0) != 0)  
        {  
            mpz_mul(mult,primes[i],primes[i]);
            mpz_set(testBoucle,mult);
            mpz_set(j,mult);
            while(mpz_cmp_ui(j,k) <0) 
            {  
                
                mpz_sub_ui(minus,j,1);
                int index=mpz_get_ui(minus);
      
                mpz_set_ui(primes[index],0);
                mpz_add(j,j,primes[i]);
            }  
        }  
        ++i;
        mpz_mul(mult,primes[i],primes[i]);
        mpz_set(testBoucle,mult);
    }  
  

}



void fact_p_1_pollard(mpz_t N,uint16_t b1,uint64_t b2)
{

 
    // first primality test  in case N is prime
    if (mpz_probab_prime_p(N, 10) == 0){
            // N is already prime 
    } else{
            printf("************* p-1 pollard ***************** \n");
            factor f;
            f.prime_factors = NULL;
            f.exp = NULL;
            f.nb_factors = 0;

            /*mpz_t fact1,fact2,B1,B2;
            double duration;
            clock_t start, finish;
            mpz_inits(fact1,fact2,B1,B2,NULL);
            mpz_set_ui(B1,b1);
            mpz_set_ui(B2,b2);
            start = clock(); // start clock 
            int r = p_minus_1(N,fact1,B1,B2);
            finish = clock(); // end clock 

            // end of algorithm printout the result 
            mpz_divexact (fact2, N, fact1);
            gmp_printf(" %Zu = %Zu  * %Zu\n",N,fact1,fact2);
            duration = (double)(finish - start)/CLOCKS_PER_SEC;
            printf("Running time: %f seconds.\n", duration);
            printf("CPU time %ld\n",cputime());*/

    }     
    
   

    

}

int fact_trialDivision(mpz_t *tab,mpz_t n,mpz_t pmax)
{
    mpz_t r,p,quotient,result;
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
	if (i!=0){ return 1;} else{	return -1;}

}
