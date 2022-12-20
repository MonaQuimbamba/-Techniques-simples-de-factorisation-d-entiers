#include "Includes/factorization.h"




int main(int argc, char const *argv[])
{


    if( argc < 2 )
   {
          printf("Usage : %s  mode   \n", argv[0]);
          printf(" mode <1> :  Divisions successives \n");
          printf(" mode <2> :  rho de Pollard  \n");
          printf(" mode <3> :  p-1  de Pollard  \n"); 
          exit(-1);
   }


   int mode = atoi(argv[1]);
   mpz_t n; 
   int r;
   factor f;
   clock_t start_time, end_time;

   switch (mode)
   {
         case 1 :
               printf("Divisions successives \n");
               mpz_t p_max;
               mpz_inits(n,p_max,NULL);
      
               f.prime_factors = NULL;
               f.exp = NULL;
               f.nb_factors = 0;
               // input n, p
               printf("Enter n: ");
               gmp_scanf("%Zd", n);
               printf("Enter p_max: ");
               gmp_scanf("%Zd", p_max);
               start_time = clock();
               r= fact_trialDivision(n,p_max,&f);
               end_time = clock();
               if(r==0)    printf("Factorization completed \n");
               else if(r==-1)  printf("Factorization incompleted  try again with a p_max  bigger \n");
               print_primes_factors(&f,n);
               printf("Process finished in %.9f secs \n", (double)(end_time-start_time)/CLOCKS_PER_SEC);
               mpz_clears(n,p_max, NULL);
               
         break;

         case 2:
            printf("rho de Pollard  \n");
            
         break;
         case 3:
           printf(" p- 1 de Pollard  \n");
           mpz_t B1,B2;
           mpz_inits(n,B1,B2,NULL);
           f.prime_factors = NULL;
           f.exp = NULL;
           f.nb_factors = 0;
            // input n, p,B1,B2
            printf("Enter n: ");
            gmp_scanf("%Zd", n);
            printf("Enter B1: ");
            gmp_scanf("%Zd", B1);
            printf("Enter B2: ");
            gmp_scanf("%Zd", B2);
            start_time = clock();
            r = fact_p_1_pollard(n,B1,B2,&f);
            end_time = clock();
            if(r==0)    printf("Factorization completed \n");
            else if(r==-1)  printf("Factorization incompleted  try again with a bound B2 bigger \n");
            print_primes_factors(&f,n);
            printf("Process finished in %.9f secs \n", (double)(end_time-start_time)/CLOCKS_PER_SEC);



            mpz_clears(n, B1, B2, NULL);
         break;


      }
 
    
   

    
    

    return 0;
}
