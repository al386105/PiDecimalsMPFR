#ifndef CHUDNOVSKY_V2_MPI
#define CHUDNOVSKY_V2_MPI

void Chudnovsky_algorithm_v2_MPI(int num_procs, int proc_id, mpfr_t pi, 
                                int num_iterations, int num_threads, int precision_bits);

#endif

