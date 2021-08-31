#include <stdio.h>
#include <stdlib.h>
#include "../../Headers/Sequential/PiCalculator.h"
#include "../../Headers/Common/Print_title.h"


int incorrect_params(char* exec_name){
    printf("  Number of params are not correct. Try with:\n");
    printf("    %s algorithm precision \n", exec_name);
    printf("\n");
}

int main(int argc, char **argv){    

    print_PiDecimals_title();
    printf("  Sequential version! \n");
    printf("\n");

    //Check the number of parameters are correct
    if(argc != 3){
        incorrect_params(argv[0]);
        exit(-1);
    }

    //Take algorithm and precision from params
    int algorithm = atoi(argv[1]);    
    int precision = atoi(argv[2]);

    calculate_Pi(algorithm, precision);

    exit(0);
}
