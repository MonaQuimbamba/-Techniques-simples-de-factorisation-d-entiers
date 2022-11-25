
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gmp.h"
#include "Includes/p_1_pollard.h"
#include "Includes/rho_pollard.h"
#include "Includes/trial_division.h"


int main(int argc, char const *argv[])
{
    /* code */
    hellop_1_pollard();
    hellorho_pollard();
    hellortrial_division();
    return 0;
}
