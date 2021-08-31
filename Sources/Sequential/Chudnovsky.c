#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>
#include <omp.h>

#define A 13591409
#define B 545140134
#define C 640320
#define D 426880
#define E 10005


/************************************************************************************
 * Miguel Pardo Navarro. 17/07/2021                                                 *
 * Chudnovsky formula implementation                                                *
 * This version does not computes all the factorials                                *
 * It implements a single-threaded method and another that can use multiple threads *
 *                                                                                  *
 ************************************************************************************
 * Chudnovsky formula:                                                              *
 *     426880 sqrt(10005)                 (6n)! (545140134n + 13591409)             *
 *    --------------------  = SUMMATORY( ----------------------------- ),  n >=0    *
 *            pi                            (n!)^3 (3n)! (-640320)^3n               *
 *                                                                                  *
 * Some operands of the formula are coded as:                                       *
 *      dep_a_dividend = (6n)!                                                      *
 *      dep_a_divisor  = (n!)^3 (3n)!                                               *
 *      e              = 426880 sqrt(10005)                                         *
 *                                                                                  *
 ************************************************************************************
 * Chudnovsky formula dependencies:                                                 *
 *                     (6n)!         (12n + 10)(12n + 6)(12n + 2)                   *
 *      dep_a(n) = --------------- = ---------------------------- * dep_a(n-1)      *
 *                 ((n!)^3 (3n)!)              (n + 1)^3                            *
 *                                                                                  *
 *      dep_b(n) = (-640320)^3n = (-640320)^3(n-1) * (-640320)^3)                   *
 *                                                                                  *
 *      dep_c(n) = (545140134n + 13591409) = dep_c(n - 1) + 545140134               *
 *                                                                                  *
 ************************************************************************************/


/*
 * An iteration of Chudnovsky formula
 */
void Chudnovsky_iteration(mpfr_t pi, int n, mpfr_t dep_a, mpfr_t dep_b, 
                            mpfr_t dep_c, mpfr_t aux){
    mpfr_mul(aux, dep_a, dep_c, MPFR_RNDN);
    mpfr_div(aux, aux, dep_b, MPFR_RNDN);
    
    mpfr_add(pi, pi, aux, MPFR_RNDN);
}

/*
 * Sequential Pi number calculation using the Chudnovsky algorithm
 * Single thread implementation
 */
void Chudnovsky_algorithm(mpfr_t pi, int num_iterations){
    int i, factor_a;
    mpfr_t dep_a, dep_a_dividend, dep_a_divisor, dep_b, dep_c, e, c, aux;
    
    struct timeval t1, t2;
    double execution_time;
    
    mpfr_inits(dep_a_dividend, dep_a_divisor, aux, NULL);
    mpfr_inits(dep_a, dep_b, dep_c, e, c, NULL);
    mpfr_set_ui(dep_a, 1, MPFR_RNDN);
    mpfr_set_ui(dep_b, 1, MPFR_RNDN);
    mpfr_set_ui(dep_c, A, MPFR_RNDN);
    mpfr_set_ui(e, E, MPFR_RNDN);
    mpfr_set_ui(c, C, MPFR_RNDN);
    mpfr_neg(c, c, MPFR_RNDN);
    mpfr_pow_ui(c, c, 3, MPFR_RNDN);

    for(i = 0; i < num_iterations; i ++){
        gettimeofday(&t1, NULL);
        Chudnovsky_iteration(pi, i, dep_a, dep_b, dep_c, aux);
        //Update dep_a:
        factor_a = (12 * i);
        mpfr_set_ui(dep_a_dividend, factor_a + 2, MPFR_RNDN); //ESTO ES UN SOBRECOSTE BRUTAL (ES POR EL + n)???!!

        mpfr_mul_ui(dep_a_dividend, dep_a_dividend, factor_a + 6, MPFR_RNDN);
        mpfr_mul_ui(dep_a_dividend, dep_a_dividend, factor_a + 10, MPFR_RNDN);
        mpfr_mul(dep_a_dividend, dep_a_dividend, dep_a, MPFR_RNDN);

        mpfr_set_ui(dep_a_divisor, i + 1, MPFR_RNDN); //ESTO ES UN SOBRECOSTE BRUTAL!!
        mpfr_pow_ui(dep_a_divisor, dep_a_divisor , 3, MPFR_RNDN);
        mpfr_div(dep_a, dep_a_dividend, dep_a_divisor, MPFR_RNDN);

        //Update dep_b:
        mpfr_mul(dep_b, dep_b, c, MPFR_RNDN);

        //Update dep_c:
        mpfr_add_ui(dep_c, dep_c, B, MPFR_RNDN);
        
        gettimeofday(&t2, NULL);
        execution_time = ((t2.tv_sec - t1.tv_sec) * 1000000u +  t2.tv_usec - t1.tv_usec)/1.e3; 
        //printf(" %f \n", execution_time);

    }
    mpfr_sqrt(e, e, MPFR_RNDN);
    mpfr_mul_ui(e, e, D, MPFR_RNDN);
    mpfr_div(pi, e, pi, MPFR_RNDN);    
    
    //Clear memory
    mpfr_clears(dep_a, dep_a_dividend, dep_a_divisor, dep_b, dep_c, e, c, aux, NULL);

}
