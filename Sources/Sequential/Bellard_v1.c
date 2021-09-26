#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include <omp.h>


/************************************************************************************
 * Miguel Pardo Navarro. 17/07/2021                                                 *
 * Bellard formula implementation                                                   *
 * It implements a single-threaded method and another that can use multiple threads *
 *                                                                                  *
 ************************************************************************************
 * Bellard formula:                                                                 *
 *                 (-1)^n     32     1      256     64       4       4       1      *
 * 2^6 * pi = SUM( ------ [- ---- - ---- + ----- - ----- - ----- - ----- + -----])  *
 *                 1024^n    4n+1   4n+3   10n+1   10n+3   10n+5   10n+7   10n+9    *
 *                                                                                  *
 * Formula quotients are coded as:                                                  *
 *             32          1           256          64                              *
 *        a = ----,   b = ----,   c = -----,   d = -----,                           *
 *            4n+1        4n+3        10n+1        10n+3                            *
 *                                                                                  *
 *              4            4            1         (-1)^n                          *
 *        e = -----,   f = -----,   g = -----,   m = -----,                         *
 *            10n+5        10n+7        10n+9        2^10n                          *
 *                                                                                  *
 ************************************************************************************
 * Bellard formula dependencies:                                                    *
 *                           1            1                                         *
 *              dep_m(n) = ------ = -----------------                               *
 *                         1024^n   1024^(n-1) * 1024                               *
 *                                                                                  *
 *              dep_a(n) = 4n  = dep_a(n-1) + 4                                     *
 *                                                                                  *
 *              dep_b(n) = 10n = dep_a(n-1) + 10                                    *
 *                                                                                  *
 ************************************************************************************/


/*
 * An iteration of Bellard formula
 */
void Bellard_iteration_v1(mpfr_t pi, int n, mpfr_t m, mpfr_t a, mpfr_t b, mpfr_t c, mpfr_t d, 
                    mpfr_t e, mpfr_t f, mpfr_t g, mpfr_t aux, int dep_a, int dep_b){
    mpfr_set_ui(a, 32, MPFR_RNDN);              // a = ( 32 / ( 4n + 1))
    mpfr_set_ui(b, 1, MPFR_RNDN);               // b = (  1 / ( 4n + 3))
    mpfr_set_ui(c, 256, MPFR_RNDN);             // c = (256 / (10n + 1))
    mpfr_set_ui(d, 64, MPFR_RNDN);              // d = ( 64 / (10n + 3))
    mpfr_set_ui(e, 4, MPFR_RNDN);               // e = (  4 / (10n + 5))
    mpfr_set_ui(f, 4, MPFR_RNDN);               // f = (  4 / (10n + 7))
    mpfr_set_ui(g, 1, MPFR_RNDN);               // g = (  1 / (10n + 9))
    mpfr_set_ui(aux, 0, MPFR_RNDN);             // aux = (- a - b + c - d - e - f + g)  

    mpfr_div_ui(a, a, dep_a + 1, MPFR_RNDN);    // a = ( 32 / ( 4n + 1))
    mpfr_div_ui(b, b, dep_a + 3, MPFR_RNDN);    // b = (  1 / ( 4n + 3))

    mpfr_div_ui(c, c, dep_b + 1, MPFR_RNDN);    // c = (256 / (10n + 1))
    mpfr_div_ui(d, d, dep_b + 3, MPFR_RNDN);    // d = ( 64 / (10n + 3))
    mpfr_div_ui(e, e, dep_b + 5, MPFR_RNDN);    // e = (  4 / (10n + 5))
    mpfr_div_ui(f, f, dep_b + 7, MPFR_RNDN);    // f = (  4 / (10n + 7))
    mpfr_div_ui(g, g, dep_b + 9, MPFR_RNDN);    // g = (  1 / (10n + 9))

    // aux = (- a - b + c - d - e - f + g)   
    mpfr_neg(a, a, MPFR_RNDN);
    mpfr_sub(aux, a, b, MPFR_RNDN);
    mpfr_sub(c, c, d, MPFR_RNDN);
    mpfr_sub(c, c, e, MPFR_RNDN);
    mpfr_sub(c, c, f, MPFR_RNDN);
    mpfr_add(c, c, g, MPFR_RNDN);
    mpfr_add(aux, aux, c, MPFR_RNDN);

    // aux = m * aux
    mpfr_mul(aux, aux, m, MPFR_RNDN);   

    mpfr_add(pi, pi, aux, MPFR_RNDN); 
}

/*
 * Sequential Pi number calculation using the Bellard algorithm
 * Single thread implementation
 */
void Bellard_algorithm_v1(mpfr_t pi, int num_iterations){   
    int i, dep_a, dep_b;
    mpfr_t dep_m, jump, a, b, c, d, e, f, g, aux;    

    dep_a = 0, dep_b = 0;       
    mpfr_init_set_d(jump, 1, MPFR_RNDN);            // jump = 1/1024  
    mpfr_div_ui(jump, jump, 1024, MPFR_RNDN); 
    mpfr_init_set_ui(dep_m, 1, MPFR_RNDN);          // dep_m = ((-1)^n)/1024)
    mpfr_inits(a, b, c, d, e, f, g, aux, NULL);

    for(i = 0; i < num_iterations; i++){ 
        Bellard_iteration_v1(pi, i, dep_m, a, b, c, d, e, f, g, aux, dep_a, dep_b);   
        // Update dependencies for next iteration: 
        mpfr_mul(dep_m, dep_m, jump, MPFR_RNDN);
        mpfr_neg(dep_m, dep_m, MPFR_RNDN);
        dep_a += 4;
        dep_b += 10;
    }

    mpfr_div_ui(pi, pi, 64, MPFR_RNDN);
    
    mpfr_clears(dep_m, jump, a, b, c, d, e, f, g, aux, NULL);
}

