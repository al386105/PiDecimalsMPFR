#!/bin/bash

program="$1"
error=""
RED_OUTPUT="tput setaf 1"
GREEN_OUTPUT="tput setaf 2"
RESET_OUTPUT="tput sgr0"

errors(){
    echo "params are not correct. They should be: ./compile.sh program"
    echo "  if program is Sequential -> compile sequential version of PiDecimalsMPFR"
    echo "  if program is OMP -> compile parallel OMP version of PiDecimalsMPFR "
    echo "  if program is MPI -> compile parallel bybrid OMP and MPI version of PiDecimalsMPFR "
    exit 1
}

#CHECK PARAMS
if [ "$#" -ne 1 ]; then
   errors
fi

if [ "$program" = "Sequential" ]; then
	error=$(gcc -o sequential.x Sources/Sequential/*.c Sources/Common/*.c -lmpfr -lgmp 2>&1 1>/dev/null)

elif [ "$program" = "OMP" ]; then
	error=$(gcc -fopenmp -o parallelOMP.x Sources/OMP/*.c Sources/Sequential/BBP*.c Sources/Sequential/Bellard*.c Sources/Sequential/Chudnovsky*.c Sources/Common/*.c -lmpfr -lgmp 2>&1 1>/dev/null)

elif [ "$program" = "MPI" ]; then 
	error=$(mpicc -fopenmp -o parallelMPI.x Sources/MPI/*.c Sources/Sequential/BBP*.c Sources/Sequential/Bellard*.c Sources/Sequential/Chudnovsky*.c Sources/Common/*.c -lmpfr -lgmp -lm 2>&1 1>/dev/null)
else
    errors
fi

#GIVE FEEDBACK ABOUT COMPILATION
if [[ -z "$error" ]]; then
    echo -n "COMPILATION "
    ${GREEN_OUTPUT}
    echo -n "DONE"
    ${RESET_OUTPUT}
    echo " SUCCESFULLY "

else
    ${RED_OUTPUT}
    echo -n "ERROR: "
    ${RESET_OUTPUT}
    echo "$error"
fi

