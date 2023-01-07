#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>

int algo_rho_de_pollard_opt(mpz_t n, mpz_t pmax, mpz_t d)
{
	 mpz_t a, b,tmp,result;
	mpz_inits(a,b,result,tmp,NULL);
	mpz_set_ui(a,2);
	mpz_set_ui(b,2);
	int i=0;
	printf("prime=%d\n",mpz_probab_prime_p(n,10) );
	if (mpz_probab_prime_p(n,10)==0)
	{
		do
		{
			mpz_set(tmp,a);

			do
			{		
				mpz_powm_ui(a,a,2,n);
				mpz_add(a,a,b);
				mpz_powm_ui(b,b,2,n);
				mpz_powm_ui(b,b,2,n);
				mpz_sub(result,a,b);
				mpz_mod(result,result,n);
				mpz_gcd(d,result,n);
				// gmp_printf("d=%Zd\n",d);
				if (mpz_cmp_ui(d,1)>0 && mpz_cmp(d,n)!=0 )
				{
					return 1;
								}
				i=i+1;
				printf("i=%d\n",i );
			}while(mpz_cmp_ui(pmax,i)>=0);
			i=0;
			mpz_add_ui(a,tmp,1);
			mpz_add_ui(b,tmp,1);
			gmp_printf("a=%Zd",a);

		
		}while(mpz_cmp_ui(d,1)<=0 || mpz_cmp(d,n)==0 );

		return 1;

	}
	else
	{
		printf("L'entier est premier:\n");
		return -1;
	}
	
}	


int main(int argc, char const *argv[])
{
	mpz_t a,b,n,d,pmax;

	mpz_inits(a,b,n,d,pmax,NULL);
	mpz_set_ui(a,2);
	mpz_set_ui(b,2);
	mpz_set_str(n,"455459293873664234764378236674635826467387352447",10);
	mpz_set_str(pmax,"50000",10);
	int res=0;
	res=algo_rho_de_pollard_opt(n,pmax,d);
	printf("result=%d\n",res );
	gmp_printf("d=%Zd",d);

	
}