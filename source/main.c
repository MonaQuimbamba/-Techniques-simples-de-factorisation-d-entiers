#include "Includes/factorization.h"



int testUnitaire(PrimeFactors *f, mpz_t n){
     int r=0;
     mpz_t test,tmp;
     mpz_inits(test,tmp,NULL);
     mpz_set_ui(test,1);
     for (int i = 0; i < f->num_factors; i++){
                 mpz_set_ui(tmp,1);
                if (f->exponents[i] != 1) {  
                    mpz_pow_ui(tmp,f->prime_factors[i],f->exponents[i]);
                }               
                else { 
                    mpz_set(tmp,f->prime_factors[i]);
                }  
                mpz_mul(test,test,tmp);
}

//gmp_printf(" n= %Zu et n computed %Zu \n ",n,test);
if(mpz_cmp(test,n)!=0) r=-1;
mpz_clears(test,tmp,NULL);

return r;

}


int main(int argc, char const *argv[])
{


    if( argc < 2 )
   {
          printf("Usage : %s  mode   \n", argv[0]);
          printf(" mode <0>  :  use all methods  with defaults parameters or instead use like this ./main <n> p_max>  <B1>  <B2>  to give your own parameters\n");
          printf(" mode <1>  :  Divisions successives \n");
          printf(" mode <2>  :  rho de Pollard  \n");
          printf(" mode <3>  :  p-1  de Pollard  \n"); 
          exit(-1);
   }

   if(argv[1]!=NULL){

        // declaration of variables 
            mpz_t n,p_max,B1,B2;
            uint64_t nb_iterations; 
            clock_t start_time, end_time;
            int r,restest;
             PrimeFactors factors;

         if( strcmp(argv[1], "1") == 0 || strcmp(argv[1], "2") == 0 || strcmp(argv[1], "3") == 0) { // switch case mode 

               int mode = atoi(argv[1]);
             
              

               switch (mode){
                     case 1 :
                           printf("Divisions successives \n");
                           mpz_inits(n,p_max,NULL);
                  
                          
                           factors.prime_factors = NULL;
                           factors.exponents = NULL;
                           factors.num_factors = 0;
                           // input n, p
                           printf("Enter n: ");
                           gmp_scanf("%Zd", n);
                           printf("Enter p_max: ");
                           gmp_scanf("%Zd", p_max);
                           start_time = clock();
                           r= fact_trialDivision(n,p_max,&factors);
                           end_time = clock();

                           restest = testUnitaire(&factors,n);
                           if(restest!=0) printf(" P -1 Pollard failled  try again with a p_max bigger \n");
                           else{

                           if(r==0)   {    printf("Factorization completed \n\n"); print_primes_factors(&factors,n); }
                           else if(r==-1)  printf("Factorization incompleted  try again with p_max bigger \n");
                           }


                           printf("Process finished in %.9f secs \n", (double)(end_time-start_time)/CLOCKS_PER_SEC);

                           // free factors  after all 
                           for (int i = 0; i < factors.num_factors; i++) mpz_clear(factors.prime_factors[i]);
                           free(factors.prime_factors);
                           free(factors.exponents); 
                           mpz_clears(n,p_max, NULL);
                           
                     break;

                     case 2:
                        printf("rho de Pollard  \n");
                        
                         mpz_inits(n,nb_iterations,NULL);
                  
                          
                           factors.prime_factors = NULL;
                           factors.exponents = NULL;
                           factors.num_factors = 0;
                           // input n, p
                           printf("Enter n: ");
                           gmp_scanf("%Zd", n);
                           printf("Enter nb iterations: ");
                           scanf("%ld", &nb_iterations);
                           start_time = clock();
                           r= fact_pollard_rho(n,&factors,nb_iterations);
                           end_time = clock();

                           restest = testUnitaire(&factors,n);
                           if(restest!=0) printf(" Rho de  Pollard failled  try again with a nb interations  bigger \n");
                           else{

                           if(r==0)   {    printf("Factorization completed \n\n"); print_primes_factors(&factors,n); }
                           else if(r==-1)  printf("Factorization incompleted  try again with nb iterations  bigger \n");
                           }


                           printf("Process finished in %.9f secs \n", (double)(end_time-start_time)/CLOCKS_PER_SEC);

                           // free factors  after all 
                           for (int i = 0; i < factors.num_factors; i++) mpz_clear(factors.prime_factors[i]);
                           free(factors.prime_factors);
                           free(factors.exponents); 
                           mpz_clears(n,p_max, NULL);
                        
                     break;
                     case 3:
                        printf(" p- 1 de Pollard  \n");
                     
                        mpz_inits(n,B1,B2,NULL);
                        factors.prime_factors = NULL;
                        factors.exponents = NULL;
                        factors.num_factors = 0;
                        // input n, p,B1,B2
                        printf("Enter n: ");
                        gmp_scanf("%Zd", n);
                        printf("Enter B1: ");
                        gmp_scanf("%Zd", B1);
                        printf("Enter B2: ");
                        gmp_scanf("%Zd", B2);
                        start_time = clock();
                        r = fact_p_1_pollard(n,B1,B2,&factors);
                        end_time = clock();

                        restest = testUnitaire(&factors,n);
                        if(restest!=0) printf(" P -1 Pollard failled  try again with a bound B1 or B2 bigger \n");
                        else{

                           if(r==0)   {    printf("Factorization completed \n\n"); print_primes_factors(&factors,n); }
                           else if(r==-1)  printf("Factorization incompleted  try again with a bound B2 bigger \n");
                        }
                        printf("Process finished in %.9f secs \n", (double)(end_time-start_time)/CLOCKS_PER_SEC);
                        for (int i = 0; i < factors.num_factors; i++) mpz_clear(factors.prime_factors[i]);
                        free(factors.prime_factors);
                        free(factors.exponents); 
                        mpz_clears(n, B1, B2, NULL);
                     break;


                  }

         }
         else if( strcmp(argv[1],"0") ==0 )
         
         {
            printf(" Running all methods  with defaults paremeters \n");
            printf("\n");
            printf(" p_max = 10^5 \n");
            printf(" B1    = 10^5 \n");
            printf(" B2    = 10^6 \n\n");

          


            // initiation of variables 
            mpz_inits(n,p_max,B1,B2,NULL);

            printf("Enter n: ");
            gmp_scanf("%Zd", n);
            //mpz_set_str(n,argv[2],10); // get N from the agrv 
            mpz_set_ui(p_max,100000); // set default value 
            mpz_set_ui(B1,100000);
            mpz_set_ui(B2,10000000);
          

          
            factors.prime_factors = NULL;
            factors.exponents = NULL;
            factors.num_factors = 0;
  

            printf("start  running Divisions successives ... \n");
            //printf("***********************************************\n");
            start_time = clock();
            r= fact_trialDivision(n,p_max,&factors);
            end_time = clock();
            // check out if the factorization was a succes
            restest = testUnitaire(&factors,n);

            if(restest!=0) printf("Divisions successives failled \n");
            else{
                     if(r==0){ printf("Factorization completed \n\n"); print_primes_factors(&factors,n);   }
                  else if(r==-1)  printf("Factorization incompleted  try again with a p_max  bigger \n");
            }
               // free factors  after all 
            for (int i = 0; i < factors.num_factors; i++) mpz_clear(factors.prime_factors[i]);
            free(factors.prime_factors);
            free(factors.exponents);         
            printf("Process finished in %.9f secs \n", (double)(end_time-start_time)/CLOCKS_PER_SEC);         
            printf("end  Divisions successives  \n\n\n");

         printf("***********************************************\n\n\n");  
     

            printf("start  running rho de Pollard   ... \n");
          
            /*factors.prime_factors = NULL;
         factors.exponents = NULL;
         factors.num_factors=0;

            // set up the defaults vars values here 
            start_time = clock();
            // run the function here 
            end_time = clock();
            // free factors  after all 
            for (int i = 0; i < f.nb_factors; i++) mpz_clear(f.prime_factors[i]);
            free(f.prime_factors);
            free(f.exp); 
            printf("Process finished in %.9f secs \n", (double)(end_time-start_time)/CLOCKS_PER_SEC);
            
            printf("end  rho de Pollard  \n\n\n");
           

            */
       

         printf("***********************************************\n");

         printf("start  p -1 de  Pollard  \n");
         // set up the defaults vars values here 
         factors.prime_factors = NULL;
         factors.exponents = NULL;
         factors.num_factors=0;
         start_time = clock();
         // run the function here 
         r = fact_p_1_pollard(n,B1,B2,&factors);
         end_time = clock();

         restest = testUnitaire(&factors,n);
         if(restest!=0) printf(" P -1 Pollard failled  try again with a bound B1 or B2 bigger \n");
         else{

         if(r==0)   {    printf("Factorization completed \n\n"); print_primes_factors(&factors,n); }
         else if(r==-1)  printf("Factorization incompleted  try again with a bound B2 bigger \n");

         }
         // free factors  after all 
         for (int i = 0; i < factors.num_factors; i++) mpz_clear(factors.prime_factors[i]);
         free(factors.prime_factors);
         free(factors.exponents); 
         printf("Process finished in %.9f secs \n", (double)(end_time-start_time)/CLOCKS_PER_SEC);

         printf("end  p -1 de  Pollard  \n");


         // clear up 
         mpz_clears(n,p_max,B1,B2,NULL);
          
            // clears up all vars 
            
         }

         else{

               // all methods with user parameters 
               printf(" news paras");

         }


   }




    
    

    return 0;
}
