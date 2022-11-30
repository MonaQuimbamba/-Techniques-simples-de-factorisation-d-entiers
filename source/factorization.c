
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
void p_1_pollard_v1(mpz_t N,mpz_t output,uint64_t b)
{
    
    mpz_t q,power_primes,a,fact,tmp_primes;
    mpz_t *primes;
    mpz_inits(q,power_primes,a,fact,tmp_primes,NULL);
    mpz_set_ui(q,1);
    mpz_set_ui(a,2);

  while(true)
  {

    
    primes = malloc(sizeof(mpz_t) * b-1);
    criblesimple(primes,b);
    //for(int i =1 ; i < b; i++) gmp_printf(" %Zu ",primes[i]);
    int i=1;
    while(i < b) {
        if(mpz_cmp_ui(primes[i],0)!=0){          
            powerPrime(primes[i],q,b);
        }
        i++;
    }
    // gmp_printf("q=%Zu b=%d\n",q,b);
    mpz_powm(power_primes,a,q,N);
    mpz_sub_ui(power_primes,power_primes,1);
    mpz_mod(power_primes,power_primes,N);
    mpz_powm_ui(fact,power_primes,1,N);
    if(mpz_cmp_ui(fact,0)!=0){
        mpz_gcd(fact,power_primes,N);
        if(mpz_cmp_ui(fact,1)==0){
            // next b 
           mpz_set_ui(tmp_primes,b);
           mpz_nextprime(tmp_primes,tmp_primes);
           b= mpz_get_ui(tmp_primes);
           mpz_set_ui(q,1);
          
        }else{
                
                if(mpz_probab_prime_p(fact,10)>0){
                    mpz_set(output,fact);
                    break;
                }
        }
      
    }

    
 }


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

void fact_p_1_pollard(mpz_t N)
{
    printf("************* p-1 pollard ***************** \n");
    mpz_t fact1,fact2;
    uint64_t b=7;
    double duration;
    clock_t start, finish;
    mpz_inits(fact1,fact2,NULL);
    printf("**********************************\n");
    printf("\n******  	version 1:		 *****\n");
    printf("\n************************************* \n");
    start = clock();
    p_1_pollard(N,fact1);
    finish = clock();
    mpz_divexact (fact2, N, fact1);
    gmp_printf("q=%Zu  p=%Zu\n",fact1,fact2);
    duration = (double)(finish - start)/CLOCKS_PER_SEC;
    printf("Running time: %f seconds.\n", duration);
    printf("CPU time %ld\n",cputime());

    printf("**********************************\n");
    printf("\n******  	version 2:		 *****\n");
    printf("\n************************************* \n");
  
    start = clock();
    p_1_pollard_v1(N,fact1,b);
    finish = clock();
    mpz_divexact (fact2, N, fact1);
    gmp_printf("q=%Zu  p=%Zu\n",fact1,fact2);
    duration = (double)(finish - start)/CLOCKS_PER_SEC;
    printf("Running time: %f seconds.\n", duration);
    printf("CPU time %ld\n",cputime());

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
		i=i-1; // on ne retourne pas le dernier facteur de n qui est superieur à pmax
	}
	if (i!=0){ return 1;} else{	return -1;}

}
