#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include <omp.h>
#include <math.h>
#include "mpi.h"
#include "../../Headers/Sequential/Bellard_v1.h"
#include "../../Headers/MPI/OperationsMPI.h"


/*
 * Parallel Pi number calculation using the Bellard algorithm
 * The number of iterations is divided by blocks, 
 * so each process calculates a part of pi using threads. 
 * Each process will cyclically divide the iterations 
 * among the threads to calculate its part.  
 * Finally, a collective reduction operation will be performed
 * using a user defined function in OperationsMPI. 
 */
void Bellard_algorithm_MPI(int num_procs, int proc_id, mpfr_t pi, 
                                int num_iterations, int num_threads, int precision_bits){
    int block_size, block_start, block_end, position, packet_size, d_elements;
    mpfr_t local_proc_pi, ONE;

    block_size = (num_iterations + num_procs - 1) / num_procs;
    block_start = proc_id * block_size;
    block_end = block_start + block_size;
    if (block_end > num_iterations) block_end = num_iterations;

    mpfr_inits2(precision_bits, ONE, local_proc_pi, NULL);
    mpfr_set_ui(local_proc_pi, 0, MPFR_RNDN);
    mpfr_set_ui(ONE, 1, MPFR_RNDN); 

    //Set the number of threads 
    omp_set_num_threads(num_threads);

    #pragma omp parallel 
    {
        int thread_id, i, dep_a, dep_b, jump_dep_a, jump_dep_b, next_i;
        mpfr_t local_thread_pi, dep_m, a, b, c, d, e, f, g, aux;

        thread_id = omp_get_thread_num();

        mpfr_init2(local_thread_pi, precision_bits);               // private thread pi
        mpfr_set_ui(local_thread_pi, 0, MPFR_RNDN);
        dep_a = (block_start + thread_id) * 4;
        dep_b = (block_start + thread_id) * 10;
        jump_dep_a = 4 * num_threads;
        jump_dep_b = 10 * num_threads;
        mpfr_init2(dep_m, precision_bits);
        mpfr_mul_2exp(dep_m, ONE, 10 * (block_start + thread_id), MPFR_RNDN);
        mpfr_div(dep_m, ONE, dep_m, MPFR_RNDN);
        if((thread_id + block_start) % 2 != 0) mpfr_neg(dep_m, dep_m, MPFR_RNDN);                 
        mpfr_inits2(precision_bits, a, b, c, d, e, f, g, aux, NULL);

        //First Phase -> Working on a local variable
        #pragma omp parallel for 
            for(i = block_start + thread_id; i < block_end; i+=num_threads){
                Bellard_iteration_v1(local_thread_pi, i, dep_m, a, b, c, d, e, f, g, aux, dep_a, dep_b);
                // Update dependencies for next iteration:
                next_i = i + num_threads;
                mpfr_mul_2exp(dep_m, ONE, 10 * next_i, MPFR_RNDN);
                mpfr_div(dep_m, ONE, dep_m, MPFR_RNDN);
                if (next_i % 2 != 0) mpfr_neg(dep_m, dep_m, MPFR_RNDN);
                dep_a += jump_dep_a;
                dep_b += jump_dep_b;  
            }

        //Second Phase -> Accumulate the result in the global variable
        #pragma omp critical
        mpfr_add(local_proc_pi, local_proc_pi, local_thread_pi, MPFR_RNDN);

        //Clear thread memory
        mpfr_free_cache();
        mpfr_clears(local_thread_pi, dep_m, a, b, c, d, e, f, g, aux, NULL);   
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
        mpfr_div_ui(pi, pi, 64, MPFR_RNDN);
    }

    //Clear memory
    MPI_Op_free(&add_op);
    mpfr_clears(local_proc_pi, ONE, NULL);       

}

