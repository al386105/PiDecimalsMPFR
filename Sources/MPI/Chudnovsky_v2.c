#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include <omp.h>
#include <math.h>
#include "mpi.h"
#include "../../Headers/Sequential/Chudnovsky_v2.h"
#include "../../Headers/MPI/OperationsMPI.h"

#define A 13591409
#define B 545140134
#define C 640320
#define D 426880
#define E 10005


/*
 * This method is used by Chudnovsky_algorithm_MPI threads
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
 * Parallel Pi number calculation using the Chudnovsky algorithm
 * The number of iterations is divided by blocks 
 * so each process calculates a part of pi with multiple threads (or just one thread). 
 * Each process will also divide the iterations in blocks
 * among the threads to calculate its part.  
 * Finally, a collective reduction operation will be performed 
 * using a user defined function in OperationsMPI. 
 */
void Chudnovsky_algorithm_v2_MPI(int num_procs, int proc_id, mpfr_t pi, 
                                    int num_iterations, int num_threads, int precision_bits){
    int block_size, block_start, block_end, position, packet_size, d_elements;
    mpfr_t local_proc_pi, e, c;

    block_size = (num_iterations + num_procs - 1) / num_procs;
    block_start = proc_id * block_size;
    block_end = block_start + block_size;
    if (block_end > num_iterations) block_end = num_iterations;

    mpfr_inits2(precision_bits, local_proc_pi, e, c, NULL);
    mpfr_set_ui(local_proc_pi, 0, MPFR_RNDN);
    mpfr_set_ui(e, E, MPFR_RNDN); 
    mpfr_set_ui(c, C, MPFR_RNDN); 
    mpfr_neg(c, c, MPFR_RNDN);
    mpfr_pow_ui(c, c, 3, MPFR_RNDN);

    //Set the number of threads 
    omp_set_num_threads(num_threads);

    #pragma omp parallel 
    {
        int thread_id, i, thread_block_size, thread_block_start, thread_block_end, factor_a;
        mpfr_t local_thread_pi, dep_a, dep_a_dividend, dep_a_divisor, dep_b, dep_c, aux;

        thread_id = omp_get_thread_num();
        thread_block_size = (block_size + num_threads - 1) / num_threads;
        thread_block_start = (thread_id * thread_block_size) + block_start;
        thread_block_end = thread_block_start + thread_block_size;
        if (thread_block_end > block_end) thread_block_end = block_end;

        mpfr_init2(local_thread_pi, precision_bits);    // private thread pi
        mpfr_set_ui(local_thread_pi, 0, MPFR_RNDN);
        mpfr_inits2(precision_bits, dep_a, dep_b, dep_c, dep_a_dividend, dep_a_divisor, aux, NULL);
        init_dep_a(dep_a, thread_block_start, precision_bits);
        mpfr_pow_ui(dep_b, c, thread_block_start, MPFR_RNDN);
        mpfr_set_ui(dep_c, B, MPFR_RNDN);
        mpfr_mul_ui(dep_c, dep_c, thread_block_start, MPFR_RNDN);
        mpfr_add_ui(dep_c, dep_c, A, MPFR_RNDN);
        factor_a = 12 * thread_block_start;


        //First Phase -> Working on a local variable        
        #pragma omp parallel for 
            for(i = thread_block_start; i < thread_block_end; i++){
                Chudnovsky_iteration(local_thread_pi, i, dep_a, dep_b, dep_c, aux);
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
        mpfr_add(local_proc_pi, local_proc_pi, local_thread_pi, MPFR_RNDN);

        //Clear thread memory
        mpfr_free_cache();
        mpfr_clears(local_thread_pi, dep_a, dep_a_dividend, dep_a_divisor, dep_b, dep_c, aux, NULL); 
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
        mpfr_sqrt(e, e, MPFR_RNDN);
        mpfr_mul_ui(e, e, D, MPFR_RNDN);
        mpfr_div(pi, e, pi, MPFR_RNDN); 
    }

    //Clear memory
    MPI_Op_free(&add_op);
    mpfr_clears(local_proc_pi, e, c, NULL);       

}