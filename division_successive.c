#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>

int division_successive(mpz_t *tab,mpz_t n,mpz_t pmax)
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
	if (i!=0)
	{
		return 1;
	}
	else
	{
		return -1;
	}

}

int main(int argc, char const *argv[])
{

	int result=0;
	mpz_t n,pmax;
	mpz_init_set_ui(n,12345);
	mpz_init_set_ui(pmax,13);
	mpz_t *tab;
	result=division_successive(tab,n,pmax);
	printf("result= %d\n",result );
	

}