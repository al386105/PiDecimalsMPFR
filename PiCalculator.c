#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include "BBPAlgorithm.h"
#include "BellardAlgorithm.h"
#include "ChudnovskyAlgorithm.h"


void checkDecimals(mpfr_t pi){
    //Cast the number we want to check to string
   
    int bytes_of_pi = ((pi -> _mpfr_prec) * sizeof(mp_limb_t));
    char calculated_pi[bytes_of_pi]; 
    mpfr_sprintf(calculated_pi, "%Re", pi);

    //Read the correct pi number from numeroPiCorrecto.txt file and compares the decimals to calculated pi
    FILE * file;
    file = fopen("./resources/numeroPiCorrecto.txt", "r");
    if(file == NULL){
        printf("numeroPiCorrecto.txt not found \n");
        exit(-1);
    } 

    char correct_pi_char;
    int i = 0;
    while((correct_pi_char = fgetc(file)) != EOF){
        if( (i >= bytes_of_pi) || (correct_pi_char != calculated_pi[i])){
            break;
        }
        i++;
    }
    i = (i < 2) ? 0: i - 2;
    printf("Match the first %d decimal places \n", i);
    
    fclose(file);

    //mpfr_printf("Pi: %Re \n", pi);
}

void checkMemoryError(mpfr_t pi){
    if (pi == NULL){
        printf("ERROR EN LA RESERVA DE MEMORIA. PRECISION DEMASIADO GRANDE \n");
        exit(0);
    }
}

void BBPAlgorithm(int num_threads, int precision){
    int num_iterations, precision_bits;
    mpfr_t pi;

    precision_bits = precision * 8;
    num_iterations = precision * 0.84;

    //Set gmp float precision (in bits) and init pi
    mpfr_set_default_prec(precision_bits); 
    mpfr_init_set_ui(pi, 0, MPFR_RNDN);

    if(num_threads <= 1){ 
        SequentialBBPAlgorithm(pi, num_iterations);
    } else {
        ParallelBBPAlgorithm(pi, num_iterations, num_threads, precision_bits);
    }
    
    checkDecimals(pi);
    
    mpfr_clear(pi);
}

void BellardAlgorithm(int num_threads, int precision){
    int num_iterations, precision_bits;
    mpfr_t pi;
    
    precision_bits = precision * 8;
    num_iterations = precision / 3;

    //Set gmp float precision (in bits) and init pi
    mpfr_set_default_prec(precision_bits); 
    mpfr_init_set_ui(pi, 0, MPFR_RNDN);

    if(num_threads <= 1){ 
        SequentialBellardAlgorithm(pi, num_iterations);
    } else {
        ParallelBellardAlgorithm(pi, num_iterations, num_threads, precision_bits);
    }
    
    checkDecimals(pi);

    mpfr_clear(pi);
}   


void ChudnovskyAlgorithm(int num_threads, int precision){
    int num_iterations, precision_bits;
    mpfr_t pi;

    precision_bits = precision * 8;
    num_iterations = (precision + 14 - 1) / 14;  //Division por exceso

    //Set gmp float precision (in bits) and init pi
    mpfr_set_default_prec(precision_bits); 
    mpfr_init_set_ui(pi, 0, MPFR_RNDN);

    if(num_threads <= 1){ 
        SequentialChudnovskyAlgorithm(pi, num_iterations);
    } else {
        ParallelChudnovskyAlgorithm(pi, num_iterations, num_threads, precision_bits);
    }
    
    checkDecimals(pi);

    mpfr_clear(pi);
}   
