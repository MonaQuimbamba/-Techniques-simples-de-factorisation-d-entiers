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
				mpz_set(fact->tab[fact->taille],result);
				mpz_sqrt(r,n);
				e=0;
				fact->taille=fact->taille+1;
				fact->tab=realloc(fact->tab,(fact->taille+1)*sizeof(mpz_t)); // reallocation dynamique de la taille du tableau
				// fact->tab[fact->taille]=malloc(sizeof(mpz_t));
				mpz_init_set_ui(fact->tab[fact->taille],1);
				printf("t=%d\n",fact->taille );
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

int main(int argc, char const *argv[])
{
	printf("omar\n");
	facteur* p;
	mpz_t n, pmax;
	mpz_init_set_ui(n,1235);
	mpz_init_set_ui(pmax,19);
	p=algo_division_succ(n,pmax);
	printf("result=%d\n",p->result);
	printf("taille=%d\n",p->taille);
	for (int i = 0; i < p->taille; ++i)
	{
		gmp_printf("tab[%d]=%Zd\n",i,p->tab[i]);
	}
	
}