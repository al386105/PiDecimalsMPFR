#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>
#include <omp.h>

int main(int argc, char **argv){ 
    mpfr_t a, b;
    mpfr_init_set_ui(a, 5, MPFR_RNDN);
    mpfr_init_set_ui(b, 6, MPFR_RNDN);
    mpfr_add(a, a, b, MPFR_RNDN);
    mpfr_printf("a: %.17Rg \n", a);

}   
