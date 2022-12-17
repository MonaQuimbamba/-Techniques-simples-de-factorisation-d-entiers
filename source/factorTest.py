#:/bin/python3
import os
import subprocess

# Using readlines()
file1 = open('40_digits_test.txt', 'r')
lines = file1.readlines()


# Strips the newline character
for i in  range(0,len(lines)):
    cmd = subprocess.Popen("./main %s 100000 100000"%lines[i].strip(), shell=True,stdout=subprocess.PIPE)
    (resultat, ignorer) = cmd.communicate()
    print(resultat)
 