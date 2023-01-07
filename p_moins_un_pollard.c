#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
typedef struct facteur facteur;
struct facteur
{
	mpz_t* tab;
	int taille;
	int result;

};

facteur* algo_division_succ(mpz_t n,mpz_t pmax)
{
	facteur* fact=(facteur*)malloc(sizeof(facteur));
	fact->taille=0;
	fact->result=0;
	fact->tab=(mpz_t*)malloc(sizeof(mpz_t));


	mpz_t r,p,quotient,result;
	mpz_inits(quotient,r,result,NULL);
	mpz_sqrt(r,n);
	mpz_init_set_ui(p,2);
	int e=0;


	
	// fact->tab[fact->taille]=(mpz_t*)malloc(sizeof(mpz_t)); //allocation dynamique du tableau des facteurs de n
	mpz_init_set_ui(fact->tab[fact->taille],1);

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
			mpz_pow_ui(result,p,e);
			if (mpz_cmp(p,pmax)<=0) // on teste si le facteur trouvé est inferieur ou pas à pmax
			{
				mpz_set(fact->tab[fact->taille],p);
				mpz_sqrt(r,n);
				e=0;
				fact->taille=fact->taille+1;
				fact->tab=realloc(fact->tab,(fact->taille+1)*sizeof(mpz_t)); // reallocation dynamique de la taille du tableau
				// fact->tab[fact->taille]=malloc(sizeof(mpz_t));
				mpz_init_set_ui(fact->tab[fact->taille],1);
				// printf("t=%d\n",fact->taille );
			}
		}
		mpz_add_ui(p,p,1); // on augmente p de 1		
	}
	mpz_init_set_ui(fact->tab[fact->taille],1);
	if(mpz_cmp_ui(n,1)>0 && mpz_cmp(n,pmax)<=0)
	{
		mpz_set(fact->tab[fact->taille],n);
		fact->taille=fact->taille+1;
	}
	

	if (fact->taille!=0)
	{
		fact->result=1;
		return fact;
	}
	else
	{
		fact->result=-1;
		return fact;
	}

}


void algo_p_1_de_pollard(mpz_t n, mpz_t B,mpz_t d)
{
	int r=0,i=0;
	mpz_t a,pmax,A;
	mpz_inits(a,pmax,A,NULL);
	mpz_init_set_ui(pmax,10000);
	facteur* fact;
	fact=algo_division_succ(B,pmax);
	gmp_printf("taille=%d\n",fact->taille);
	r=fact->taille;
	do
	{
		
		while(mpz_cmp_ui(a,0)==0 || mpz_cmp(a,n)==0 )
		{
			mpz_inits(a,NULL);
			gmp_randstate_t alea;
			gmp_randinit_default(alea); 
			gmp_randseed_ui(alea, time(NULL));
			mpz_urandomm (a, alea,n);
			gmp_randclear (alea);
		}
		gmp_printf("a=%Zd\n",a);
		
		mpz_set(A,a);
			mpz_sub_ui(A,A,1);
			mpz_gcd(d,n,A);
		while(mpz_cmp_ui(d,1)==0 )
		{
			printf("omar____\n");
			printf("i=%d\n",i);
			mpz_set(A,a);
			mpz_sub_ui(A,A,1);
			mpz_gcd(d,n,A);
			mpz_powm(a,a,fact->tab[i],n);
			gmp_printf("d=%Zd\n",d);
			i=i+1;
			if(i==fact->taille)
			{
				break;
			}
		}
		i=0;
		printf("okkkkkkkkkkkkkkkkkkkk\n");

	}while(mpz_cmp(d,n)==0 || mpz_cmp_ui(d,1)==0);


}


int main(int argc, char const *argv[])
{
	mpz_t n,B,d;
	mpz_init(d);
	mpz_init_set_ui(n,1123456789);
	mpz_init_set_ui(B,100);
	algo_p_1_de_pollard(n, B,d);
	gmp_printf("d final=%Zd",d);
	
}