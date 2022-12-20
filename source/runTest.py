#:/bin/python3
import os
import subprocess
import argparse
import numpy as np


"""# Using readlines()
file1 = open('40_digits_test.txt', 'r')
lines = file1.readlines()


# Strips the newline character
for i in  range(0,len(lines)):
    cmd = subprocess.Popen("./main %s 100000 100000"%lines[i].strip(), shell=True,stdout=subprocess.PIPE)
    (resultat, ignorer) = cmd.communicate()
    print(resultat)
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

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
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
        

