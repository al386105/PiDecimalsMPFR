#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>
#include <omp.h>
#include "../../Headers/Sequential/Chudnovsky.h"


#define A 13591409
#define B 545140134
#define C 640320
#define D 426880
#define E 10005


/*
 * This method is used by ParallelChudnovskyAlgorithm threads
 * for computing the first value of dep_a
 */
void init_dep_a(mpfr_t dep_a, int block_start, int precision_bits){
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
 * This method provides an optimal distribution for each thread
 * based on the Chudnovsky iterations analysis.
 * It returns an array of three integers:
 *   distribution[0] -> block size
 *   distribution[1] -> block start
 *   distribution[2] -> block end 
 */
int * get_thread_distribution(int num_threads, int thread_id, int num_iterations){
    int * distribution, i, block_size, block_start, block_end, row, column; 
    FILE * ratios_file;
    float working_ratios[160][41];

    //Open the working_ratios file 
    ratios_file = fopen("Resources/working_ratios.txt", "r");
    if(ratios_file == NULL){
        printf("working_ratios.txt not found \n");
        exit(-1);
    } 

    //Load the working_ratios matrix 
    row = 0;
    while (fscanf(ratios_file, "%f", &working_ratios[row][0]) == 1){
        for (column = 1; column < 41; column++){
            fscanf(ratios_file, "%f", &working_ratios[row][column]);
        }
        row++;
    }

    distribution = malloc(sizeof(int) * 3);
    if(num_threads == 1){
        distribution[0] = num_iterations;
        distribution[1] = 0;
        distribution[2] = num_iterations;
        return distribution; 
    }

    block_size = working_ratios[thread_id][num_threads / 4] * num_iterations / 100;
    block_start = 0;
    for(i = 0; i < thread_id; i ++){
        block_start += working_ratios[i][num_threads / 4] * num_iterations / 100;
    }
    block_end = block_start + block_size;
    if (thread_id == num_threads -1) block_end = num_iterations;
    
    distribution[0] = block_size;
    distribution[1] = block_start;
    distribution[2] = block_end;

    return distribution;
}


/*
 * Parallel Pi number calculation using the Chudnovsky algorithm
 * Multiple threads can be used
 * The number of iterations is divided by blocks 
 * so each thread calculates a part of pi.  
 */
void Chudnovsky_algorithm_OMP(mpfr_t pi, int num_iterations, int num_threads, int precision_bits){
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
        int thread_id, i, block_size, block_start, block_end, factor_a, * distribution;
        mpfr_t local_pi, dep_a, dep_a_dividend, dep_a_divisor, dep_b, dep_c, aux;

        thread_id = omp_get_thread_num();
//        distribution = get_thread_distribution(num_threads, thread_id, num_iterations);
//        block_size = distribution[0];
//        block_start = distribution[1];
//        block_end = distribution[2];
//
        
        block_size = (num_iterations + num_threads - 1) / num_threads;
        block_start = thread_id * block_size;
        block_end = block_start + block_size;
        if (block_end > num_iterations) block_end = num_iterations;
        
        mpfr_inits2(precision_bits, local_pi, dep_a, dep_b, dep_c, dep_a_dividend, dep_a_divisor, aux, NULL);
        mpfr_set_ui(local_pi, 0, MPFR_RNDN);    // private thread pi
        init_dep_a(dep_a, block_start, precision_bits);
        mpfr_pow_ui(dep_b, c, block_start, MPFR_RNDN);
        mpfr_set_ui(dep_c, B, MPFR_RNDN);
        mpfr_mul_ui(dep_c, dep_c, block_start, MPFR_RNDN);
        mpfr_add_ui(dep_c, dep_c, A, MPFR_RNDN);
        factor_a = 12 * block_start;

        //First Phase -> Working on a local variable        
        #pragma omp parallel for 
            for(i = block_start; i < block_end; i++){
                Chudnovsky_iteration(local_pi, i, dep_a, dep_b, dep_c, aux);
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