#include <stdio.h>
#include <gmp.h>

// Finds all prime factors of the composite number n within the specified bounds using the Pollard-p-1 algorithm
void pollard_p_1(mpz_t n, mpz_t lower_bound, mpz_t upper_bound) {
    // Initialize the variables p, B, and g
    mpz_t p, B, g;
    mpz_inits(p, B, g, NULL);

    // Set p to the lower bound
    mpz_set(p, lower_bound);

    // Loop until the upper bound is reached
    while (mpz_cmp(p, upper_bound) <= 0) {
        // Calculate B = p^(e+1) + 1, where e is the largest exponent such that B is not divisible by n
        int e = 1;
        mpz_pow_ui(B, p, e+1);
        mpz_add_ui(B, B, 1);
        while (mpz_divisible_p(n, B) == 0) {
            e++;
            mpz_pow_ui(B, p, e+1);
            mpz_add_ui(B, B, 1);
        }

        // Calculate g = gcd(B-1, n)
        mpz_sub_ui(B, B, 1);
        mpz_gcd(g, B, n);

        // If g is not equal to 1 or n, then g is a factor of n
        if (mpz_cmp_ui(g, 1) != 0 && mpz_cmp(g, n) != 0) {
            gmp_printf("Found factor: %Zd\n", g);
            // Divide n by g and set n to the result
            mpz_divexact(n, n, g);
        }

        // Increment p and try again
        mpz_add_ui(p, p, 1);
    }

    // Clear the variables
    mpz_clears(p, B, g, NULL);
}

int main() {
    // Test the Pollard-p-1 algorithm with a composite number and bounds
    mpz_t n, lower_bound, upper_bound;
    mpz_inits(n, lower_bound, upper_bound, NULL);
    mpz_set_str(n, "15", 10);
    mpz_set_ui(lower_bound, 2);
    mpz_set_ui(upper_bound, 5);
    pollard_p_1(n, lower_bound, upper_bound);
    mpz_clears(n, lower_bound, upper_bound, NULL);
    return 0;
}