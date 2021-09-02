#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include <time.h>
#include "../../Headers/OMP/BBP.h"
#include "../../Headers/OMP/Bellard.h"
#include "../../Headers/OMP/Chudnovsky.h"
#include "../../Headers/Common/Check_decimals.h"


double gettimeofday();

void check_errors_OMP(int precision, int num_iterations, int num_threads){
    if (precision <= 0){
        printf("  Precision should be greater than cero. \n\n");
        exit(-1);
    } 
    if (num_iterations < num_threads){
        printf("  The number of iterations required for the computation is too small to be solved with %d threads. \n", num_threads);
        printf("  Try using a greater precision or lower threads number. \n\n");
        exit(-1);
    }
}

void print_running_properties_OMP(int precision, int num_iterations, int num_threads){
    printf("  Precision used: %d \n", precision);
    printf("  Iterations done: %d \n", num_iterations);
    printf("  Number of threads: %d\n", num_threads);
}

void calculate_Pi_OMP(int algorithm, int precision, int num_threads){
    double execution_time;
    struct timeval t1, t2;
    mpfr_t pi;
    int num_iterations, decimals_computed, precision_bits;

    precision_bits = 8 * precision;
    
    gettimeofday(&t1, NULL);

    //Set mpfr float precision (in bits) and init pi
    mpfr_set_default_prec(precision_bits); 
    mpfr_init_set_ui(pi, 0, MPFR_RNDN);
    
    switch (algorithm)
    {
    case 0:
        num_iterations = precision * 0.84;
        check_errors_OMP(precision, num_iterations, num_threads);
        printf("  Algorithm: BBP \n");
        print_running_properties_OMP(precision, num_iterations, num_threads);
        BBP_algorithm_OMP(pi, num_iterations, num_threads, precision_bits);
        break;

    case 1:
        num_iterations = precision / 3;
        check_errors_OMP(precision, num_iterations, num_threads);
        printf("  Algorithm: Bellard \n");
        print_running_properties_OMP(precision, num_iterations, num_threads);
        Bellard_algorithm_OMP(pi, num_iterations, num_threads, precision_bits);
        break;
    
    case 2:
        num_iterations = (precision + 14 - 1) / 14;  //Division por exceso
        check_errors_OMP(precision, num_iterations, num_threads);
        printf("  Algorithm: Chudnovsky (Last version) \n");
        print_running_properties_OMP(precision, num_iterations, num_threads);
        Chudnovsky_algorithm_OMP(pi, num_iterations, num_threads, precision_bits);
        break;
    
    default:
        printf("  Algorithm selected is not correct. Try with: \n");
        printf("      algorithm == 0 -> BBP  \n");
        printf("      algorithm == 1 -> Bellard \n");
        printf("      algorithm == 2 -> Chudnovsky  \n");
        printf("\n");
        exit(-1);
        break;
    }

    gettimeofday(&t2, NULL);
    execution_time = ((t2.tv_sec - t1.tv_sec) * 1000000u +  t2.tv_usec - t1.tv_usec)/1.e6; 
    decimals_computed = check_decimals(pi);
    mpfr_clear(pi);
    printf("  Match the first %d decimals \n", decimals_computed);
    printf("  Execution time: %f seconds \n", execution_time);
    printf("\n");
}

