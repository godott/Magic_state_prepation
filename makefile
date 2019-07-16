CC=mpicc
FFLAGS=-I.
OMPFLAGS=-fopenmp
CCFLAGS=-std=gnu99
CFLAGS=$(OMPFLAGS) $(FFLAGS) $(CCFLAGS)
DEPS=error_propagation.h parser.h
OBJ=search_flag.o error_propagation.o parser.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

flagsearch: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)
