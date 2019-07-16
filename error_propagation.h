#ifndef ERROR_PROPAGATION_H
#define ERROR_PROPAGATION_H

#include <stdio.h>
#include <stdlib.h>
#include <regex.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>

/* Macros */

#define MAX_FILENAME_SIZE 100
#define MAX_QUBIT_NUM 70
#define MAX_GATENAME_LEN 6
#define MAX_GATE_OPS 2
#define MAX_INSTRUCTION_ON_QUBIT 100
#define MAX_INSTRUCTION_NUM 800
#define INVALID_SIZE -1
#define MAX_ERROR_NUM 12

/* Data structure and type */
typedef struct {
	char name[MAX_GATENAME_LEN];
	size_t op[MAX_GATE_OPS];
	size_t qcirc_ind[MAX_GATE_OPS];
} instruction;  // For each gate, it has a name, qubits to operate on and it has an id on each qubit.

/* Error data type. A: abnormal; N: not specified. */
typedef enum error_t {I, X, Y, Z, IX, IY, IZ, 
                      XI, XX, XY, XZ, 
                      YI, YX, YY, YZ,
                      ZI, ZX, ZY, ZZ, 
                      H, HX, HY, HZ, 
                      A, N} error_t;

typedef struct {
	size_t len[MAX_QUBIT_NUM];
    error_t result[MAX_QUBIT_NUM][MAX_ERROR_NUM];
} raw_error; // error array, with number of errors

typedef struct {
	instruction ins_id[MAX_INSTRUCTION_NUM];
    size_t qcirc[MAX_QUBIT_NUM][MAX_INSTRUCTION_ON_QUBIT];
    size_t len_qcirc[MAX_QUBIT_NUM];
    size_t currentQubitNum;
    size_t dataQubitNum;
    size_t flagQubitNum;
    size_t ancillaInd;     
    size_t len_ins;
} flag_circ; // Quantum circuit structure


typedef struct{
    error_t syndrome[MAX_QUBIT_NUM];
    size_t loc[15];
    error_t error_type[15];
    size_t loc_len;
} syndrome_ht;

typedef struct{
    syndrome_ht synHT[MAX_INSTRUCTION_NUM * 3];
    size_t len;
} syndrome_ind;

/* function declaration */
void clifford_propagation(error_t error_pauli, size_t qbit_index, size_t depth, 
                          flag_circ *flagCirc, raw_error * rawError);

void error_gen(size_t gate_id, error_t error_type, flag_circ * flagCirc, raw_error * rawError);

void print_raw_error_result(size_t qubitNum, raw_error * rawError);

void error_to_string(error_t raw, char error_string[]);

void raw_to_final_error_result(size_t qubitNum, raw_error * rawError, error_t error_result[]);

error_t combine_error(error_t first, error_t second);
#endif
