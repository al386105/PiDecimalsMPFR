#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "../../Headers/MPI/PiCalculator.h"
#include "../../Headers/Common/Print_title.h"


int incorrect_params(char* exec_name){
    printf("  Number of params are not correct. Try with:\n");
    printf("    mpirun -np num_procs %s algorithm precision num_threads\n", exec_name);
    printf("\n");
}

int main(int argc, char **argv){    
    int num_procs, proc_id;

    //Init MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id); 

    //Print PiDecimals title 
    if(proc_id == 0){
        print_PiDecimals_title();
        printf("  This version is done for clusters!\n");
        printf("\n");
    }

    //Check the number of parameters are correct
    if(argc != 4){
        incorrect_params(argv[0]);
        exit(-1);
    }

    //Take operation, precision and number of threads from params
    int algorithm = atoi(argv[1]);    
    int precision = atoi(argv[2]);
    int num_threads = (atoi(argv[3]) <= 0) ? 1 : atoi(argv[3]);

    //Compute Pi
    calculate_Pi_MPI(num_procs, proc_id, algorithm, precision, num_threads);

    MPI_Finalize();

    exit(0);
}