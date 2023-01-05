/*
    Author: Nhat BUI
    Email : van-nhat.bui@etu.unilim.fr
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "time.h"
#include "gmp.h"
#include "math.h"
#include <string.h>

#define testNum 50000000

int pollard_rho_Floy_cycle(mpz_t n, mpz_t d, gmp_randstate_t my_generator){

    mpz_t t,x, y, c;
    mpz_inits(t,x, y, c, NULL);

    // Set the initial values for x and y
    mpz_set_ui(x, 2);
    mpz_set_ui(y, 2);

    mpz_urandomm(c,my_generator,n);
    

    // Set the initial value for d
    mpz_set_ui(d, 1);
    unsigned int i = 0;
    while(mpz_cmp_ui(d, 1) == 0){
        if (i > testNum) break;
        i++;

        // Set x to (x^2 + c) mod n
        mpz_mul(x, x, x);
        mpz_add(x, x, c);
        mpz_mod(x, x, n);

        // Set y to (y^2 + c) mod n
        mpz_mul(y, y, y);
        mpz_add(y, y, c);
        mpz_mod(y, y, n);

        // Set y to (y^2 + c) mod n again
        mpz_mul(y, y, y);
        mpz_add(y, y, c);
        mpz_mod(y, y, n);
        
        // d = gcd(|x-y|,n)
        mpz_sub(d, x, y);
        mpz_abs(d, d);
        mpz_gcd(d, d, n);

        if ( mpz_cmp(d,n)!=0 && mpz_cmp_ui(d,1) !=0) {
            break;
        }
    }
    
    // Print the result
    if (mpz_cmp(d,n)==0 || mpz_cmp_ui(d,1)==0 ){
        gmp_printf("Can't factorize number %Zd\n",n);

    } else{
        gmp_printf("One of the divisors for %Zd is %Zd\n",n,d);

        mpz_divexact(t,n,d);
        mpz_mod(x,n,d);
        gmp_printf("%Zd = %Zd * %Zd + %Zd\n",n,d,t,x);
    }
    mpz_clears(t,x, y, c, NULL);
}      

int pollard_rho_Brent_cycle(mpz_t n, mpz_t d, gmp_randstate_t my_generator){
    mpz_t x,y,c,t;
    mpz_inits(x,c,y,t,NULL);
    long long int r, k, m;

    // set random for y, c, m
    mpz_urandomm(c,my_generator,n);

    // initialize y and p
    mpz_set_ui(y, 2);
    mpz_set_ui(d, 1);

    // r is a power of 2
    r = 1;
    
    while ((mpz_cmp_ui(d, 1) == 0) || (mpz_cmp(d, n) == 0)){
        
        // x to the current position of y, namely x_i
        mpz_set(x, y);

        k = 0;

        // y to x_(i+r)
        for (int i=0;i++;i<r){
            mpz_mul(y, y, y);
            mpz_add(y, y, c);
            mpz_mod(y, y, n);
        }

        // test gcd(|x-y|, n); where y ranges from x_(i+r) to x_(i+2r)
        while ((k < r) && (mpz_cmp_ui(d, 1) == 0)){
            mpz_mul(y, y, y);
            mpz_add(y, y, c);
            mpz_mod(y, y, n);
            mpz_sub(t, x, y);
            mpz_abs(t, t);
            mpz_mod(t, t, n);                  
            mpz_gcd(d, t, n);
            k++;
        }

        // going to the next power of 2
        r = r*2;
        if (r >= testNum) break;
        if ( mpz_cmp(d,n)!=0 && mpz_cmp_ui(d,1) !=0) {
            // printf("r = %lld\n",r);
            break;
        }
    }

    // Print the result
    if (mpz_cmp(d,n)==0 || mpz_cmp_ui(d,1)==0 ){
        gmp_printf("Can't factorize number %Zd\n",n);

    } else{
        gmp_printf("One of the divisors for %Zd is %Zd\n",n,d);

        mpz_divexact(t,n,d);
        mpz_mod(x,n,d);
        gmp_printf("%Zd = %Zd * %Zd + %Zd\n",n,d,t,x);
    }
    
    mpz_clears(x,c,y,t,NULL);

    return 1;
}

void run(mpz_t n, mpz_t d, gmp_randstate_t my_generator){
    clock_t timeBegin, timeEnd;
    double timeSpent;
    // FILE* time_measurement;
    
    gmp_printf("Factorize number %Zd \n", n);
    printf("------------------------------------\n");

    printf("Run by Pollard rho with Floy cycle detection.\n");
    timeBegin = clock();
    pollard_rho_Floy_cycle(n,d,my_generator);
    timeEnd = clock();
    timeSpent = (float)(timeEnd - timeBegin) / CLOCKS_PER_SEC;
    printf("Executed time: %.4f\n", timeSpent);

    printf("------------------------------------\n");

    printf("Run by Pollard rho with Brent cycle detection.\n");
    timeBegin = clock();
    pollard_rho_Brent_cycle(n,d,my_generator);
    timeEnd = clock();
    timeSpent = (float)(timeEnd - timeBegin) / CLOCKS_PER_SEC;
    printf("Executed time: %.4f\n", timeSpent);

}

int main(int argc, char* argv[]){

    char* input1 = "1125939825397831601";
    char* input2 = "925276410789441750962080530947";
    char* input3 = "115792089237316195423570985008687907853269984665640564039457584007913129639937";

    // Initialize variables
    mpz_t n,d;
    mpz_inits(n,d,NULL);

    // random state
    gmp_randstate_t my_generator; 
    gmp_randinit_default(my_generator);
    gmp_randseed_ui(my_generator, time(NULL));

    if (argc != 2)
    {
        printf("Usage : %s Number\n", argv[0]);
        printf("If no input, run default test case.\n");
        printf("\n");

        printf("************* Test 1 ***************\n");
        mpz_set_str(n, input1,10);
        run(n,d,my_generator);
        printf("************************************\n");
        printf("\n");

        printf("************* Test 2 ***************\n");
        mpz_set_str(n, input2,10);
        run(n,d,my_generator);
        printf("************************************\n");
        printf("\n");

        printf("************* Test 3 ***************\n");
        mpz_set_str(n, input3,10);
        run(n,d,my_generator);
        printf("************************************\n");

    } else {
        // Get input
        printf("************* Run ***************\n");
        mpz_set_str(n, argv[1],10);
        run(n,d,my_generator);
        printf("************************************\n");

    }

    gmp_randclear(my_generator);
    mpz_clears(n,d,NULL);
}