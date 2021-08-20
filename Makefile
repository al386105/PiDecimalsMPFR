GCC = gcc

LIBS = -lmpfr -lgmp 

EXECS = parallel.x

default: all

all: $(EXECS)

parallel.x :  PiDecimals.o PiCalculator.o ChudnovskyAlgorithm.o BellardAlgorithm.o BBPAlgorithm.o 
	$(GCC) -fopenmp -o parallel.x PiDecimals.o PiCalculator.o ChudnovskyAlgorithm.o BellardAlgorithm.o BBPAlgorithm.o $(LIBS)

.c.o:
	$(GCC) -fopenmp -c $*.c

clean:
	rm *.o

clear:
	rm *.o $(EXECS)

# gcc -fopenmp -o parallel.x PiDecimals.c PiCalculator.c ChudnovskyAlgorithm.c BellardAlgorithm.c BBPAlgorithm.c lmpfr -lgmp