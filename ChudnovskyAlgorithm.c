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
void ChudnovskyIteration(mpfr_t pi, int n, mpfr_t dep_a, mpfr_t dep_b, 
                            mpfr_t dep_c, mpfr_t aux){
    mpfr_mul(aux, dep_a, dep_c, MPFR_RNDN);
    mpfr_div(aux, aux, dep_b, MPFR_RNDN);
    
    mpfr_add(pi, pi, aux, MPFR_RNDN);
}

/*
 * Sequential Pi number calculation using the Chudnovsky algorithm
 * Single thread implementation
 */
void SequentialChudnovskyAlgorithm(mpfr_t pi, int num_iterations){
    int i, factor_a;
    mpfr_t dep_a, dep_a_dividend, dep_a_divisor, dep_b, dep_c, e, c, aux;
    
    struct timeval t1, t2;
    double execution_time;
    
    mpfr_inits2(240000, dep_a_dividend, dep_a_divisor, aux, NULL);
    mpfr_inits2(240000, dep_a, dep_b, dep_c, e, c, NULL);
    mpfr_set_ui(dep_a, 1, MPFR_RNDN);
    mpfr_set_ui(dep_b, 1, MPFR_RNDN);
    mpfr_set_ui(dep_c, A, MPFR_RNDN);
    mpfr_set_ui(e, E, MPFR_RNDN);
    mpfr_set_ui(c, C, MPFR_RNDN);
    mpfr_neg(c, c, MPFR_RNDN);
    mpfr_pow_ui(c, c, 3, MPFR_RNDN);

    for(i = 0; i < num_iterations; i ++){
        gettimeofday(&t1, NULL);
        ChudnovskyIteration(pi, i, dep_a, dep_b, dep_c, aux);
        //Update dep_a:
        factor_a = (12 * i);
        mpfr_set_ui(dep_a_dividend, factor_a + 2, MPFR_RNDN); //ESTO ES UN SOBRECOSTE BRUTAL (ES POR EL + 10)???!!

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
    printf("Iterations done: %d\n", i);
    mpfr_sqrt(e, e, MPFR_RNDN);
    mpfr_mul_ui(e, e, D, MPFR_RNDN);
    mpfr_div(pi, e, pi, MPFR_RNDN);    
    
    //Clear memory
    mpfr_clears(dep_a, dep_a_dividend, dep_a_divisor, dep_b, dep_c, e, c, aux, NULL);

}

/*
 * This method is used by ParallelChudnovskyAlgorithm threads
 * for computing the first value of dep_a
 */
void initDepA(mpfr_t dep_a, int block_start, int precision_bits){
    mpz_t factorial_n, dividend, divisor;
    mpfr_t float_dividend, float_divisor;
    mpz_inits(factorial_n, dividend, divisor, NULL);
    mpfr_inits2(precision_bits ,float_dividend, float_divisor, NULL);

    mpz_fac_ui(factorial_n, block_start);
    mpz_fac_ui(divisor, 3 * block_start);
    mpz_fac_ui(dividend, 6 * block_start);

    mpz_pow_ui(factorial_n, factorial_n, 3);
    mpz_mul(divisor, divisor, factorial_n);

    mpfr_set_z(float_dividend, dividend, MPFR_RNDN);
    mpfr_set_z(float_divisor, divisor, MPFR_RNDN);

    mpfr_div(dep_a, float_dividend, float_divisor, MPFR_RNDN);

    mpz_clears(factorial_n, dividend, divisor, NULL);
    mpfr_clears(float_dividend, float_divisor, NULL);
}


/*
 * Parallel Pi number calculation using the Chudnovsky algorithm
 * Multiple threads can be used
 * The number of iterations is divided by blocks 
 * so each thread calculates a part of pi.  
 */
void ParallelChudnovskyAlgorithm(mpfr_t pi, int num_iterations, int num_threads, int precision_bits){
    mpfr_t e, c;

    mpfr_inits2(precision_bits, e, c, NULL);
    mpfr_set_ui(e, E, MPFR_RNDN);
    mpfr_set_ui(c, C, MPFR_RNDN);
    mpfr_neg(c, c, MPFR_RNDN);
    mpfr_pow_ui(c, c, 3, MPFR_RNDN);

    //Set the number of threads 
    omp_set_num_threads(num_threads);

    #pragma omp parallel 
    {   
        int thread_id, i, block_size, block_start, block_end, factor_a;
        mpfr_t local_pi, dep_a, dep_a_dividend, dep_a_divisor, dep_b, dep_c, aux;

        thread_id = omp_get_thread_num();
        block_size = (num_iterations + num_threads - 1) / num_threads;
        block_start = thread_id * block_size;
        block_end = (thread_id == num_threads - 1) ? num_iterations : block_start + block_size;
        
        mpfr_inits2(precision_bits, local_pi, dep_a, dep_b, dep_c, dep_a_dividend, dep_a_divisor, aux, NULL);
        mpfr_set_ui(local_pi, 0, MPFR_RNDN);    // private thread pi
        initDepA(dep_a, block_start, precision_bits);
        mpfr_pow_ui(dep_b, c, block_start, MPFR_RNDN);
        mpfr_set_ui(dep_c, B, MPFR_RNDN);
        mpfr_mul_ui(dep_c, dep_c, block_start, MPFR_RNDN);
        mpfr_add_ui(dep_c, dep_c, A, MPFR_RNDN);
        factor_a = 12 * block_start;

        //First Phase -> Working on a local variable        
        #pragma omp parallel for 
            for(i = block_start; i < block_end; i++){
                ChudnovskyIteration(local_pi, i, dep_a, dep_b, dep_c, aux);
                //Update dep_a:
                mpfr_set_ui(dep_a_dividend, factor_a + 10, MPFR_RNDN);
                mpfr_mul_ui(dep_a_dividend, dep_a_dividend, factor_a + 6, MPFR_RNDN);
                mpfr_mul_ui(dep_a_dividend, dep_a_dividend, factor_a + 2, MPFR_RNDN);
                mpfr_mul(dep_a_dividend, dep_a_dividend, dep_a, MPFR_RNDN);

                mpfr_set_ui(dep_a_divisor, i + 1, MPFR_RNDN);
                mpfr_pow_ui(dep_a_divisor, dep_a_divisor , 3, MPFR_RNDN);
                mpfr_div(dep_a, dep_a_dividend, dep_a_divisor, MPFR_RNDN);
                factor_a += 12;

                //Update dep_b:
                mpfr_mul(dep_b, dep_b, c, MPFR_RNDN);

                //Update dep_c:
                mpfr_add_ui(dep_c, dep_c, B, MPFR_RNDN);
            }

        //Second Phase -> Accumulate the result in the global variable 
        #pragma omp critical
        mpfr_add(pi, pi, local_pi, MPFR_RNDN);
        
        //Clear thread memory
        mpfr_clears(local_pi, dep_a, dep_b, dep_c, dep_a_dividend, dep_a_divisor, aux, NULL);   
    }

    mpfr_sqrt(e, e, MPFR_RNDN);
    mpfr_mul_ui(e, e, D, MPFR_RNDN);
    mpfr_div(pi, e, pi, MPFR_RNDN);    
    
    //Clear memory
    mpfr_clears(c, e, NULL);
}