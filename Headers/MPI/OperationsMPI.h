#ifndef OPERATIONS_MPI
#define OPERATIONS_MPI

void add(void *, void *, int *, MPI_Datatype *);
void mul(void *, void *, int *, MPI_Datatype *);
int pack(void *, mpfr_t);
void unpack(void *, mpfr_t);

#endif
