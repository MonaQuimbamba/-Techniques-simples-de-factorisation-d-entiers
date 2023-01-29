#!/bin/python3
import os
import subprocess
import argparse
import numpy as np
import random

# Using readlines()
"""file1 = open('40_digits_test.txt', 'r')
lines = file1.readlines()


# Strips the newline character
for i in  range(0,len(lines)):
    cmd = subprocess.Popen("./main %s 100000 100000"%lines[i].strip(), shell=True,stdout=subprocess.PIPE)
    (resultat, ignorer) = cmd.communicate()
    print(resultat)
    4 digits :  90866896475539
 18 digits :  569540548505160227
 20 digits :  83773924823808287865
 30 digits :  954225503633997496066854182573
 40 digits :  8831840195884534949755106936208903789520
 50 digits :  34702918589293429985076500088325701653957366943359
 """


# Generate a random float with 50 decimal places
"""random_float = np.random.uniform(0, 1)
n  = f"{random_float:.50f}"
n=n[2:]
print(len(n), " et ",n)
random_float = np.random.uniform(0, 1)
n  = f"{random_float:.40f}"
n=n[2:]
print(len(n), " et ",n)"""

def replace_zeros(numbers):
    for i in range(len(numbers)):
        if numbers[i] == 0:
            numbers[i] = random.randint(1, 100)
    return "".join(map(str, numbers))

if __name__ == "__main__":
    num_digits=60
    n = str(random.randint(10**(num_digits-1), 10**num_digits-1))
    print(" 60 digits : ",n)
    cmd = subprocess.Popen("./main %s 1000000000 100000000 1000000000 10000000000"%n, shell=True,stdout=subprocess.PIPE)
    (resultat, ignorer) = cmd.communicate()
    print(resultat)

    num_digits=65
    n = str(random.randint(10**(num_digits-1), 10**num_digits-1))
    print(" 65 digits : ",n)
    cmd = subprocess.Popen("./main %s 1000000000 100000000 1000000000 10000000000"%n, shell=True,stdout=subprocess.PIPE)
    (resultat, ignorer) = cmd.communicate()
    print(resultat)

    num_digits=70
    n = str(random.randint(10**(num_digits-1), 10**num_digits-1))
    print(" 70 digits : ",n)
    cmd = subprocess.Popen("./main %s 1000000000 100000000 1000000000 10000000000"%n, shell=True,stdout=subprocess.PIPE)
    (resultat, ignorer) = cmd.communicate()
    print(resultat)



    num_digits=75
    n = str(random.randint(10**(num_digits-1), 10**num_digits-1))
    print(" 75 digits : ",n)
    cmd = subprocess.Popen("./main %s 1000000000 100000000 1000000000 10000000000"%n, shell=True,stdout=subprocess.PIPE)
    (resultat, ignorer) = cmd.communicate()
    print(resultat)

    num_digits=80
    n = str(random.randint(10**(num_digits-1), 10**num_digits-1))
    print(" 80 digits : ",n)
    cmd = subprocess.Popen("./main %s 1000000000 100000000 1000000000 10000000000"%n, shell=True,stdout=subprocess.PIPE)
    (resultat, ignorer) = cmd.communicate()
    print(resultat)

    """parser = argparse.ArgumentParser()
    parser.add_argument("-n","--numbersize", help=" the size of the number which we want test [ex: -n 40 / -n 50]  ", required=True)
    parser.add_argument('-m','--method' , help = 'The method which we want test  ex: -m 1 for Divisions successives | -m 2 for rho de Pollard  | -m 3 p -1 Pollard ' , required=True)
    args = parser.parse_args()
    argsdict = vars(args)
    if argsdict['method']=="1":
        print("Testing the trial division method ")
        p = 100000
        print("=================================================")
        print("we are using  10^5 as p_max")
        print("if you want change it just entrer the new p_max  ")
        tmp= input("if not then just press enter :")
        if tmp:
            p=int(tmp)
        size = argsdict['numbersize']
        if size =="40":
            random_float = np.random.uniform(0, 1)
            n  = f"{random_float:.40f}"
            n=n[2:]
            cmd = subprocess.Popen("./test 1 %s %s"%(n,p), shell=True,stdout=subprocess.PIPE)
            (resultat, ignorer) = cmd.communicate()
            print(resultat)
        elif size=="50":
            random_float = np.random.uniform(0, 1)
            n  = f"{random_float:.50f}"
            n=n[2:]
            cmd = subprocess.Popen("./test 1 %s %s"%(n,p), shell=True,stdout=subprocess.PIPE)
            (resultat, ignorer) = cmd.communicate()
            print(resultat)


                    
    elif argsdict['method']=="2":
        print(" rho de Pollard ")
    elif argsdict['method']=="3":
        print(" p -1 Pollard")
        print("=========================================================")
        B1 = 100000
        B2 = 1000000
        print("=================================================")
        print("we are using  B1 = 10^5 and B2 = 10^6")
        print("if you want change B1 ?  ")
        tmp_b1= input("if not then just press enter :")
        if tmp_b1:
            B1=int(tmp_b1)

        print("if you want change B2 ?  ")
        tmp_b2= input("if not then just press enter :")
        if tmp_b2:
            B2=int(tmp_b2)

       
        size = argsdict['numbersize']
        if size =="40":
            random_float = np.random.uniform(0, 1)
            n  = f"{random_float:.40f}"
            n=n[2:]
            cmd = subprocess.Popen("./test 3 %s %s %s"%(n,B1,B2), shell=True,stdout=subprocess.PIPE)
            (resultat, ignorer) = cmd.communicate()
            print(resultat)
        elif size=="50":
            random_float = np.random.uniform(0, 1)
            n  = f"{random_float:.50f}"
            n=n[2:]
            cmd = subprocess.Popen("./test 3 %s %s %s"%(n,B1,B2), shell=True,stdout=subprocess.PIPE)
            (resultat, ignorer) = cmd.communicate()
            print(resultat)
        

"""
"""
/* *********************  TO DO ****************************************
              mpz_set(b,a); // set b(o)=a(B1)
                mpz_nextprime(p,B1);
                int i=0;
                while(mpz_cmp(B2,p)>0){ // We then compute b1 = b(o)^l1  mod n, ...


                    if(i<1){ // to compute the b1 only 
                        

                        mpz_set(tmp_b,b);// set tmp_b <- b(k-1)
                        mpz_set(tmp_p,p); // set tmp_p <- p(k-1)

                        mpz_powm(b,b,p,n); // set b <- b(k+1)
                        mpz_sub_ui(tmp, b, 1);
                        mpz_gcd(tmp, tmp, n);
                        if (mpz_cmp_ui(tmp, 1) > 0 && mpz_cmp(tmp, n) < 0){
                            mpz_set(d, tmp);
                            found = true;
                            break;
                        }
                       mpz_nextprime(p,p);
                       i++;
                    }
                    else{
                          //bk+1 = bkclk+1−lk mod n
                          mpz_sub(tmp_exp,p,tmp_p);  // compute  expo  <- p(k+1) - p(k-1)
                          mpz_powm(t,tmp_b,tmp_exp,n); // t <- b(o)^(p - tmp_p) mod n 
                          mpz_mul(t,t,b); // t <- b(k+1) * b(k-1)
                          mpz_mod(t,t,n);  // b <- t mod n    => set b <- b(k+1)

                          mpz_set(tmp_b,b); // set tmp_b <- b(k - 1)
                          mpz_set(tmp_p,p); // set tmp_p <- p(k -1)
                          mpz_set(b,t);

                        mpz_sub_ui(tmp, b, 1);
                        mpz_gcd(tmp, tmp, n);
                        if (mpz_cmp_ui(tmp, 1) > 0 && mpz_cmp(tmp, n) < 0){
                            mpz_set(d, tmp);
                            found = true;
                            break;
                        }
                       mpz_nextprime(p,p);
                    }    
                }
                mpz_clears(b,tmp_b,tmp_p,tmp_exp,t,NULL);*/
"""
