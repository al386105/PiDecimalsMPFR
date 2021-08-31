#ifndef CHUDNOVSKY_OMP
#define CHUDNOVSKY_OMP

void Chudnovsky_algorithm_OMP(mpfr_t pi, int num_iterations, int num_threads, int);
void init_dep_a(mpfr_t dep_a, int block_start, int);

#endif

