#include <stdio.h>
#include <stdlib.h>
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
 * This version computes all the factorials needed before performing the iterations *
 * It implements a single-threaded method and another that can use multiple threads *
 *                                                                                  *
 ************************************************************************************
 * Chudnovsky formula:                                                              *
 *     426880 sqrt(10005)                 (6n)! (545140134n + 13591409)             *
 *    --------------------  = SUMMATORY( ----------------------------- ),  n >=0    *
 *            pi                            (n!)^3 (3n)! (-640320)^3n               *
 *                                                                                  *
 * Some operands of the formula are coded as:                                       *
 *      dividend = (6n)! (545140134n + 13591409)                                    *
 *      divisor  = (n!)^3 (3n)! (-640320)^3n                                        *
 *      e        = 426880 sqrt(10005)                                               *
 *                                                                                  *
 ************************************************************************************
 * Chudnovsky formula dependencies:                                                 *
 *              dep_a(n) = (6n)!                                                    *
 *              dep_b(n) = (n!)^3                                                   *
 *              dep_c(n) = (3n)!                                                    *
 *              dep_d(n) = (-640320)^(3n) = (-640320)^(3 (n-1)) * (-640320)^3       *
 *              dep_e(n) = (545140134n + 13591409) = dep_c(n - 1) + 545140134       *
 *                                                                                  *
 ************************************************************************************/

/*
 * This method calculates the factorials from 0 to num_factorials (included) 
 * and stores them in their corresponding vector position (factorials[n] = n!): 
 * factorials[0] = 1, factorials[1] = 1, factorials[2] = 2, factorials[3] = 6, etc.
 * The computation is performed with a single thread. 
 */
void getFactorials(mpfr_t * factorials, int num_factorials){
    int i;
    mpfr_t f;
    mpfr_init_set_ui(f, 1, MPFR_RNDN);
    mpfr_init_set_ui(factorials[0], 1, MPFR_RNDN);
    for(i = 1; i <= num_factorials; i++){
        mpfr_mul_ui(f, f, i, MPFR_RNDN);
        mpfr_init_set(factorials[i], f, MPFR_RNDN);
    }
    mpfr_clear(f);
}

/*
 * This method clears the factorials computed and stored in mpfr_t * factorials
 */
void clearFactorials(mpfr_t * factorials, int num_factorials){
    int i;
    for(i = 0; i <= num_factorials; i++){
        mpfr_clear(factorials[i]);
    }
}

/*
 * An iteration of Chudnovsky formula
 */
void ChudnovskyIteration(mpfr_t pi, int n, mpfr_t dep_a, mpfr_t dep_b, mpfr_t dep_c, 
                        mpfr_t dep_d, mpfr_t dep_e, mpfr_t dividend, mpfr_t divisor){
    mpfr_mul(dividend, dep_a, dep_e, MPFR_RNDN);

    mpfr_mul(divisor, dep_b, dep_c, MPFR_RNDN);
    mpfr_mul(divisor, divisor, dep_d, MPFR_RNDN);
    
    mpfr_div(dividend, dividend, divisor, MPFR_RNDN);

    mpfr_add(pi, pi, dividend, MPFR_RNDN);
}

/*
 * Sequential Pi number calculation using the Chudnovsky algorithm
 * Single thread implementation
 */
void SequentialChudnovskyAlgorithm(mpfr_t pi, int num_iterations){
    int num_factorials, i; 
    num_factorials = (num_iterations * 6) + 2;
    mpfr_t factorials[num_factorials + 1];
    getFactorials(factorials, num_factorials);   

    mpfr_t dep_a, dep_b, dep_c, dep_d, dep_e, e, c, dividend, divisor;
    mpfr_inits(dividend, divisor, NULL);
    mpfr_init_set_ui(dep_a, 1, MPFR_RNDN);
    mpfr_init_set_ui(dep_b, 1, MPFR_RNDN);
    mpfr_init_set_ui(dep_c, 1, MPFR_RNDN);
    mpfr_init_set_ui(dep_d, 1, MPFR_RNDN);
    mpfr_init_set_ui(dep_e, A, MPFR_RNDN);
    mpfr_init_set_ui(e, E, MPFR_RNDN);
    mpfr_init_set_ui(c, C, MPFR_RNDN);
    mpfr_neg(c, c, MPFR_RNDN);
    mpfr_pow_ui(c, c, 3, MPFR_RNDN);

    for(i = 0; i < num_iterations; i ++){
        ChudnovskyIteration(pi, i, dep_a, dep_b, dep_c, dep_d, dep_e, dividend, divisor);
        //Update dependencies
        mpfr_set(dep_a, factorials[6 * (i + 1)], MPFR_RNDN);
        mpfr_pow_ui(dep_b, factorials[i + 1], 3, MPFR_RNDN);
        mpfr_set(dep_c, factorials[3 * (i + 1)], MPFR_RNDN);
        mpfr_mul(dep_d, dep_d, c, MPFR_RNDN);
        mpfr_add_ui(dep_e, dep_e, B, MPFR_RNDN);
    }
    mpfr_sqrt(e, e, MPFR_RNDN);
    mpfr_mul_ui(e, e, D, MPFR_RNDN);
    mpfr_div(pi, e, pi, MPFR_RNDN);    
    
    //Clear memory
    clearFactorials(factorials, num_factorials);
    mpfr_clears(dep_a, dep_b, dep_c, dep_d, dep_e, c, e, dividend, divisor, NULL);

}

/*
 * Parallel Pi number calculation using the Chudnovsky algorithm
 * Multiple threads can be used
 * The number of iterations is divided by blocks 
 * so each thread calculates a part of pi.  
 */
void ParallelChudnovskyAlgorithm(mpfr_t pi, int num_iterations, int num_threads){
    mpfr_t e, c;
    int num_factorials, block_size;
    
    num_factorials = num_iterations * 6;
    mpfr_t factorials[num_factorials + 1];
    getFactorials(factorials, num_factorials);

    block_size = (num_iterations + num_threads - 1) / num_threads;
    mpfr_init_set_ui(e, E, MPFR_RNDN);
    mpfr_init_set_ui(c, C, MPFR_RNDN);
    mpfr_neg(c, c, MPFR_RNDN);
    mpfr_pow_ui(c, c, 3, MPFR_RNDN);

    //Set the number of threads 
    omp_set_num_threads(num_threads);

    #pragma omp parallel 
    {   
        int thread_id, i, block_start, block_end;
        mpfr_t local_pi, dep_a, dep_b, dep_c, dep_d, dep_e, dividend, divisor;

        thread_id = omp_get_thread_num();
        block_start = thread_id * block_size;
        if (block_start >= num_iterations) block_start = num_iterations;
        block_end = block_start + block_size;
        if (block_end >= num_iterations) block_end = num_iterations;

        mpfr_init_set_ui(local_pi, 0, MPFR_RNDN);    // private thread pi
        mpfr_inits(dividend, divisor, NULL);
        mpfr_init_set(dep_a, factorials[block_start * 6], MPFR_RNDN);
        mpfr_init_set(dep_b, factorials[block_start], MPFR_RNDN);
        mpfr_pow_ui(dep_b, dep_b, 3, MPFR_RNDN);
        mpfr_init_set(dep_c, factorials[block_start * 3], MPFR_RNDN);
        mpfr_init_set_ui(dep_d, C, MPFR_RNDN);
        mpfr_neg(dep_d, dep_d, MPFR_RNDN);
        mpfr_pow_ui(dep_d, dep_d, block_start * 3, MPFR_RNDN);
        mpfr_init_set_ui(dep_e, B, MPFR_RNDN);
        mpfr_mul_ui(dep_e, dep_e, block_start, MPFR_RNDN);
        mpfr_add_ui(dep_e, dep_e, A, MPFR_RNDN);

        //First Phase -> Working on a local variable        
        #pragma omp parallel for 
            for(i = block_start; i < block_end; i++){
                ChudnovskyIteration(local_pi, i, dep_a, dep_b, dep_c, dep_d, dep_e, dividend, divisor);
                //Update dependencies
                mpfr_set(dep_a, factorials[6 * (i + 1)], MPFR_RNDN);
                mpfr_pow_ui(dep_b, factorials[i + 1], 3, MPFR_RNDN);
                mpfr_set(dep_c, factorials[3 * (i + 1)], MPFR_RNDN);
                mpfr_mul(dep_d, dep_d, c, MPFR_RNDN);
                mpfr_add_ui(dep_e, dep_e, B, MPFR_RNDN);
            }


        //Second Phase -> Accumulate the result in the global variable 
        #pragma omp critical
        mpfr_add(pi, pi, local_pi, MPFR_RNDN);
        
        //Clear thread memory
        mpfr_clears(local_pi, dep_a, dep_b, dep_c, dep_d, dep_e, dividend, divisor, NULL);   
    }

    mpfr_sqrt(e, e, MPFR_RNDN);
    mpfr_mul_ui(e, e, D, MPFR_RNDN);
    mpfr_div(pi, e, pi, MPFR_RNDN);    
    
    //Clear memory
    clearFactorials(factorials, num_factorials);
    mpfr_clears(c, e, NULL);
}
