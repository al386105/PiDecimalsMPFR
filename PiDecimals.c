#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "PiCalculator.h"

/*************************************************************************************************
                                                                                                
    Compiling the file: gcc -fopenmp PiCalculator.c BPPAlgorithm.c ... -o PiDecimals -lgmp       
    Executing: ./PiDecimals algorithm precision num_threads                                       
        Algorithm can be:                                                                           
        0-1 -> BBP (Bailey-Borwein-Plouffe)
        2   -> Bellard
        3-4 -> Chudnovsky
 
 *************************************************************************************************/  

double gettimeofday();

int incorrectParams(){
    printf("Introduced params are not correct. Try as: \n");
    printf(" ./PiDecimals algorithm precision num_threads \n");
    printf("    Algorithm can be: \n");
    printf("        0 -> BBP (Bailey-Borwein-Plouffe) \n");
    printf("        1 -> Bellard \n"); 
    printf("        2 -> Chudnovsky \n"); 
}

void piDecimalsTitle(){
    FILE * file;
    file = fopen("./resources/piDecimalsTitle.txt", "r");
    char character;

    int i = 0;
    while((character = fgetc(file)) != EOF){
        printf("%c", character);
    }
      
    fclose(file);
}


int main(int argc, char **argv){    
    int i;
    //Check the number of parameters are correct
    if(argc != 4){
        incorrectParams();
        exit(-1);
    }

    piDecimalsTitle();

    //Take operation, precision, number of iterations and number of threads from params
    int algorithm = atoi(argv[1]);    
    int precision = atoi(argv[2]);
    int num_threads = (atoi(argv[3]) <= 0) ? 1 : atoi(argv[3]);

    //Declare clock variables and start time
    double execution_time;
    struct timeval t1, t2; 
    gettimeofday(&t1, NULL);

    switch (algorithm)
    {
    case 0:
        printf("Algorithm: BBP \n");
        BBPAlgorithm(num_threads, precision); 
        break;
    case 1:
        printf("Algorithm: Bellard \n");
        BellardAlgorithm(num_threads, precision); 
        break;
    case 2:
        printf("Algorithm: Chudnovsky  \n");
        ChudnovskyAlgorithm(num_threads, precision); 
        break;
    default:
        incorrectParams();
        exit(-1);
        break;
    }

    
    gettimeofday(&t2, NULL);
    execution_time = ((t2.tv_sec - t1.tv_sec) * 1000000u +  t2.tv_usec - t1.tv_usec)/1.e6; 
    
    //Print the results
    printf("Precision used: %d \n", precision);
    printf("Number of threads: %d\n", num_threads);
    printf("Execution time: %f seconds. \n", execution_time);

    exit(0);
}