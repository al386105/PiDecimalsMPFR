#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include <time.h>
#include "../../Headers/Sequential/BBP.h"
#include "../../Headers/Sequential/Bellard.h"
#include "../../Headers/Sequential/Chudnovsky.h"
#include "../../Headers/Common/Check_decimals.h"


double gettimeofday();

void check_errors(int precision, int num_iterations){
    if (precision <= 0){
        printf("  Precision should be greater than cero. \n\n");
        exit(-1);
    } 
}

void print_running_properties(int precision, int num_iterations){
    printf("  Precision used: %d \n", precision);
    printf("  Iterations done: %d \n", num_iterations);
}

void calculate_Pi(int algorithm, int precision){
    double execution_time;
    struct timeval t1, t2;
    mpfr_t pi;
    int num_iterations, decimals_computed, precision_bits;
    
    precision_bits = precision * 8;
    gettimeofday(&t1, NULL);

    //Set mpfr float precision (in bits) and init pi

    mpfr_set_default_prec(precision_bits); 
    mpfr_init_set_ui(pi, 0, MPFR_RNDN);
    
    switch (algorithm)
    {
    case 0:
        num_iterations = precision * 0.84;
        check_errors(precision, num_iterations);
        printf("  Algorithm: BBP \n");
        print_running_properties(precision, num_iterations);
        BBP_algorithm(pi, num_iterations);
        break;

    case 1:
        num_iterations = precision / 3;
        check_errors(precision, num_iterations);
        printf("  Algorithm: Bellard \n");
        print_running_properties(precision, num_iterations);
        Bellard_algorithm(pi, num_iterations);
        break;
    
    case 2:
        num_iterations = (precision + 14 - 1) / 14;  //Division por exceso
        check_errors(precision, num_iterations);
        printf("  Algorithm: Chudnovsky (Last version) \n");
        print_running_properties(precision, num_iterations);
        Chudnovsky_algorithm(pi, num_iterations);
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

