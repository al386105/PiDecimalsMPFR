#ifndef BELLARD_MPI
#define BELLARD_MPI

void Bellard_algorithm_MPI(int num_procs, int proc_id, mpfr_t pi, 
                                int num_iterations, int num_threads, int precision_bits);

#endif

