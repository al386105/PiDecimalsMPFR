#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include <omp.h>
#include <math.h>
#include "mpi.h"
#include "../../Headers/Sequential/BBP.h"
#include "../../Headers/MPI/OperationsMPI.h"

#define QUOTIENT 0.0625


/*
 * Parallel Pi number calculation using the BBP algorithm
 * The number of iterations is divided by blocks, 
 * so each process calculates a part of pi using threads. 
 * Each process will also divide the iterations in blocks
 * among the threads to calculate its part.  
 * Finally, a collective reduction operation will be performed
 * using a user defined function in OperationsMPI. 
 */
void BBP_algorithm_MPI(int num_procs, int proc_id, mpfr_t pi, 
                                int num_iterations, int num_threads, int precision_bits){
    int block_size, block_start, block_end, position, packet_size, d_elements;
    mpfr_t local_proc_pi, quotient;

    block_size = (num_iterations + num_procs - 1) / num_procs;
    block_start = proc_id * block_size;
    block_end = block_start + block_size;
    if (block_end > num_iterations) block_end = num_iterations;

    mpfr_inits2(precision_bits, local_proc_pi, quotient, NULL);
    mpfr_set_d(quotient, QUOTIENT, MPFR_RNDN);
    mpfr_set_ui(local_proc_pi, 0, MPFR_RNDN);


    //Set the number of threads 
    omp_set_num_threads(num_threads);

    #pragma omp parallel 
    {
        int thread_id, i, thread_block_size, thread_block_start, thread_block_end;
        mpfr_t local_thread_pi, dep_m, quot_a, quot_b, quot_c, quot_d, aux;

        thread_id = omp_get_thread_num();
        thread_block_size = (block_size + num_threads - 1) / num_threads;
        thread_block_start = (thread_id * thread_block_size) + block_start;
        thread_block_end = thread_block_start + thread_block_size;
        if (thread_block_end > block_end) thread_block_end = block_end;
        
        mpfr_init2(local_thread_pi, precision_bits);               // private thread pi
        mpfr_set_ui(local_thread_pi, 0, MPFR_RNDN);
        mpfr_init2(dep_m, precision_bits);
        mpfr_pow_ui(dep_m, quotient, thread_block_start, MPFR_RNDN);    // m = (1/16)^n                  
        mpfr_inits2(precision_bits, quot_a, quot_b, quot_c, quot_d, aux, NULL);
        

        //First Phase -> Working on a local variable        
        #pragma omp parallel for 
            for(i = thread_block_start; i < thread_block_end; i++){
                BBP_iteration(local_thread_pi, i, dep_m, quot_a, quot_b, quot_c, quot_d, aux);
                // Update dependencies:  
                mpfr_mul(dep_m, dep_m, quotient, MPFR_RNDN);
            }

        //Second Phase -> Accumulate the result in the global variable
        #pragma omp critical
        mpfr_add(local_proc_pi, local_proc_pi, local_thread_pi, MPFR_RNDN);

        //Clear thread memory
        mpfr_free_cache();
        mpfr_clears(local_thread_pi, dep_m, quot_a, quot_b, quot_c, quot_d, aux, NULL);   
    }

    //Create user defined operation
    MPI_Op add_op;
    MPI_Op_create((MPI_User_function *)add, 0, &add_op);

    //Set buffers for cumunications and position for pack and unpack information
    d_elements = (int) ceil((float) local_proc_pi -> _mpfr_prec / (float) GMP_NUMB_BITS);
    packet_size = 8 + sizeof(mpfr_exp_t) + (d_elements * sizeof(mp_limb_t));
    char recbuffer[packet_size];
    char sendbuffer[packet_size];

    //Pack local_proc_pi in sendbuffuer
    position = pack(sendbuffer, local_proc_pi);

    //Reduce piLocal
    MPI_Reduce(sendbuffer, recbuffer, position, MPI_PACKED, add_op, 0, MPI_COMM_WORLD);

    //Unpack recbuffer in global Pi and do the last operation
    if (proc_id == 0){
        unpack(recbuffer, pi);
    }

    //Clear memory
    MPI_Op_free(&add_op);
    mpfr_clears(local_proc_pi, quotient, NULL);       

}

