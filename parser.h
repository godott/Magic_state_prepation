#ifndef PARSER_H
#define PARSER_H

#include "error_propagation.h"

/* Function neclaration */
void handle_input(int argc, char *argv[], size_t * maxE, flag_circ *flagCirc, char fileN[]);

void parse_qasm(char * fileName, flag_circ *flagCirc);
        
void parse_line(ssize_t read, char * line, flag_circ *flagCirc);

void print_summary(size_t maxError, flag_circ *flagCirc);

char *strlwr(char *str);

#endif
