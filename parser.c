#include <stdio.h>
#include <stdlib.h>
#include <regex.h>
#include <unistd.h>
#include <string.h>
#include <sys/sysinfo.h>
#include <ctype.h>
#include "parser.h"

/* The following functions are for parsing and handling input */
void handle_input(int argc, char *argv[], size_t * maxE, flag_circ *flagCirc, char fileN[]){

    int opt, fileNameLength;
    
    while((opt=getopt(argc, argv, ":h:f:e:n:")) != -1){
    
        switch(opt){
            case 'h':
                   printf("Use the following options: \
                           \n\t-n (number of flag qubits you want to try); \
                           \n\t-e (maximum number of errors, i.e., the w for w-flag circuit).\
						   \n\t-f (file name of the input .qasm)\n");
            case 'f': 
                   fileNameLength = strlen(optarg); 
                   if(fileNameLength >  MAX_FILENAME_SIZE - 1){
                        printf("File name too long! Please change the file name.");
                        exit(1);
                   }

                   strncpy(fileN, optarg, fileNameLength);
                   break;

            case 'e':
                   *maxE = atoi(optarg);
                   break;
			
			case 'n':
				   flagCirc->flagQubitNum = atoi(optarg);
				   break;

            case ':':
                   printf("The option needs a value.\n");
                   exit(1);

            case '?':
                   printf("Option not supported.\nOnly support: \
                           \n\t-n (number of flag qubits you want to try); \
                           \n\t-e (maximum number of errors).\
						   \n\t-f (file name of the input .qasm)\n");
                   exit(1);
        }
    }

	if(*maxE == INVALID_SIZE){
		printf("Please use -e to specify the maximum number of error you consider.");
		exit(1);
	}

	if(strlen(fileN) == 0){
		printf("Please use -f to specify the .qasm filename.");
		exit(1);
	}


    return;
}

void parse_qasm(char * fileName, flag_circ *flagCirc){

    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    
    fp = fopen(fileName, "r");
    if(fp == NULL){
        printf("File does not exist in current directory.\
                Please choose a valid qasm file. \n");
        exit(1);
    }

    while((read = getline(&line, &len, fp)) != -1){
        parse_line(read, line, flagCirc);
    }

    fclose(fp);
    if(line) free(line);
    
	if(flagCirc->dataQubitNum == INVALID_SIZE){
		printf("Please specify the qubit number in the .qasm file");
		exit(1);
	}

    return;
}

void parse_line(ssize_t read, char * line, flag_circ *flagCirc){

    regex_t re;
    int reti;
    size_t max_group = 4;
    regmatch_t group_array[max_group];
    
    char * reg_qubit = "qubit ([0-9]+)";
    char * reg_measz = "measz ([0-9]+)";
    char * reg_measx = "measx ([0-9]+)";
    char * reg_prepz = "prepz ([0-9]+)";
    char * reg_prepx = "prepx ([0-9]+)";
	char * reg_single = "(\\w+) ([0-9]+)";
	char * reg_two = "(\\w+) ([0-9]+),([0-9]+)";

	char sourceCopy[strlen(line) + 1];
	strcpy(sourceCopy, line);

    // Match qubit definitions
    reti = regcomp(&re, reg_qubit, REG_ICASE | REG_EXTENDED);

    if(regexec(&re, line, max_group, group_array, 0) == 0){

		size_t g = 1, dq;
		sourceCopy[group_array[g].rm_eo] = 0;
		flagCirc->dataQubitNum = atoi(sourceCopy + group_array[g].rm_so) - 1; //  read the data qubit number

		return;
    } 

	// Match measz
    reti = regcomp(&re, reg_measz, REG_ICASE | REG_EXTENDED);

    if(regexec(&re, line, max_group, group_array, 0) == 0){
		size_t g = 1, op_num;
		sourceCopy[group_array[g].rm_eo] = 0;
		op_num = atoi(sourceCopy + group_array[g].rm_so);

		strcpy(flagCirc->ins_id[flagCirc->len_ins].name, "measz");
		flagCirc->ins_id[flagCirc->len_ins].op[0] = op_num;
		flagCirc->ins_id[flagCirc->len_ins].qcirc_ind[0] = flagCirc->len_qcirc[op_num];
		
		flagCirc->qcirc[op_num][flagCirc->len_qcirc[op_num]] = flagCirc->len_ins;
		flagCirc->len_qcirc[op_num] += 1;
		flagCirc->len_ins += 1;

		return;
    } 

	// Match measx
    reti = regcomp(&re, reg_measx, REG_ICASE | REG_EXTENDED);

    if(regexec(&re, line, max_group, group_array, 0) == 0){
		size_t g = 1, op_num;
		sourceCopy[group_array[g].rm_eo] = 0;
		op_num = atoi(sourceCopy + group_array[g].rm_so);

		strcpy(flagCirc->ins_id[flagCirc->len_ins].name, "measx");
		flagCirc->ins_id[flagCirc->len_ins].op[0] = op_num;
		flagCirc->ins_id[flagCirc->len_ins].qcirc_ind[0] = flagCirc->len_qcirc[op_num];
		
		flagCirc->qcirc[op_num][flagCirc->len_qcirc[op_num]] = flagCirc->len_ins;
		flagCirc->len_qcirc[op_num] += 1;
		flagCirc->len_ins += 1;

		return;
    } 

	// Match prepz
    reti = regcomp(&re, reg_prepz, REG_ICASE | REG_EXTENDED);

    if(regexec(&re, line, max_group, group_array, 0) == 0){
		size_t g = 1, op_num;
		sourceCopy[group_array[g].rm_eo] = 0;
		op_num = atoi(sourceCopy + group_array[g].rm_so);
		
		strcpy(flagCirc->ins_id[flagCirc->len_ins].name, "prepz");
		flagCirc->ins_id[flagCirc->len_ins].op[0] = op_num;
		flagCirc->ins_id[flagCirc->len_ins].qcirc_ind[0] = flagCirc->len_qcirc[op_num];
		
		flagCirc->qcirc[op_num][flagCirc->len_qcirc[op_num]] = flagCirc->len_ins;
		flagCirc->len_qcirc[op_num] += 1;
		flagCirc->len_ins += 1;

		return;
    } 

	// Match prepx
    reti = regcomp(&re, reg_prepx, REG_ICASE | REG_EXTENDED);

    if(regexec(&re, line, max_group, group_array, 0) == 0){
		size_t g = 1, op_num;
		sourceCopy[group_array[g].rm_eo] = 0;
		op_num = atoi(sourceCopy + group_array[g].rm_so);
		
		strcpy(flagCirc->ins_id[flagCirc->len_ins].name, "prepx");
		flagCirc->ins_id[flagCirc->len_ins].op[0] = op_num;
		flagCirc->ins_id[flagCirc->len_ins].qcirc_ind[0] = flagCirc->len_qcirc[op_num];
		
		flagCirc->qcirc[op_num][flagCirc->len_qcirc[op_num]] = flagCirc->len_ins;
		flagCirc->len_qcirc[op_num] += 1;
		flagCirc->len_ins += 1;

		return;
    } 

    // Match 2-qubit gates, have to place before match 1-qubit gates
    reti = regcomp(&re, reg_two, REG_ICASE | REG_EXTENDED);

    if(regexec(&re, line, max_group, group_array, 0) == 0){
        size_t g_name = 1, g_op1 = 2, g_op2 = 3, op_num1, op_num2;
		sourceCopy[group_array[g_name].rm_eo] = 0;
		sourceCopy[group_array[g_op1].rm_eo] = 0;
		sourceCopy[group_array[g_op2].rm_eo] = 0;
		
		op_num1 = atoi(sourceCopy + group_array[g_op1].rm_so);
		op_num2 = atoi(sourceCopy + group_array[g_op2].rm_so);

		strcpy(flagCirc->ins_id[flagCirc->len_ins].name, strlwr(sourceCopy + group_array[g_name].rm_so));
		flagCirc->ins_id[flagCirc->len_ins].op[0] = op_num1;
		flagCirc->ins_id[flagCirc->len_ins].op[1] = op_num2;
		flagCirc->ins_id[flagCirc->len_ins].qcirc_ind[0] = flagCirc->len_qcirc[op_num1];
		flagCirc->ins_id[flagCirc->len_ins].qcirc_ind[1] = flagCirc->len_qcirc[op_num2];

		flagCirc->qcirc[op_num1][flagCirc->len_qcirc[op_num1]] = flagCirc->len_ins;
		flagCirc->qcirc[op_num2][flagCirc->len_qcirc[op_num2]] = flagCirc->len_ins;
		flagCirc->len_qcirc[op_num1] += 1;
		flagCirc->len_qcirc[op_num2] += 1;

		flagCirc->len_ins += 1;
		return;
    } 
	
    // Match 1 qubit gates
    reti = regcomp(&re, reg_single, REG_ICASE | REG_EXTENDED);

    if(regexec(&re, line, max_group, group_array, 0) == 0){
        size_t g_name = 1, g_op = 2, op_num;
		sourceCopy[group_array[g_name].rm_eo] = 0;
		sourceCopy[group_array[g_op].rm_eo] = 0;
		
		op_num = atoi(sourceCopy + group_array[g_op].rm_so);

		strcpy(flagCirc->ins_id[flagCirc->len_ins].name, strlwr(sourceCopy + group_array[g_name].rm_so));
		flagCirc->ins_id[flagCirc->len_ins].op[0] = op_num;
		flagCirc->ins_id[flagCirc->len_ins].qcirc_ind[0] = flagCirc->len_qcirc[op_num];
		
		flagCirc->qcirc[op_num][flagCirc->len_qcirc[op_num]] = flagCirc->len_ins;
		flagCirc->len_qcirc[op_num] += 1;

		flagCirc->len_ins += 1;
		return;
    } 

    regfree(&re);

}


char *strlwr(char *str){
  unsigned char *p = (unsigned char *)str;

  while (*p) {
     *p = tolower((unsigned char)*p);
      p++;
  }

  return str;
}


void print_summary(size_t maxError, flag_circ *flagCirc){

	printf("Summary\n");
	printf("-------------------------------\n");
	printf("-------------------------------\n");

	printf("This system has %d CPU processors configured and "
    		"%d CPU processors available.\n",
    		get_nprocs_conf(), get_nprocs());
    
	printf("-------------------------------\n");
	printf("-------------------------------\n");
    // printf("Input file name: %s\n", fileName);
    printf("Maximum error considered: %zu\n", maxError);
    printf("Flag qubit used: %zu\n", flagCirc->flagQubitNum);

	printf("data qubit number:%zu\n", flagCirc->dataQubitNum);
	printf("Total number of gates is: %zu\n", flagCirc->len_ins);

	printf("The gates are:\n");
	printf("-------------\n");

    for(size_t i = 0; i < flagCirc->len_ins; i++){
        printf("%s %zu ", flagCirc->ins_id[i].name, flagCirc->ins_id[i].op[0]);
        if(strcmp(flagCirc->ins_id[i].name, "cnot") == 0 || strcmp(flagCirc->ins_id[i].name, "ch") == 0){
            printf("%zu\n", flagCirc->ins_id[i].op[1]);
        }else{
            printf("\n");
        }
    }

	for(size_t i = 0; i < flagCirc->currentQubitNum; i++){

		printf("Number of gates(including prep and meas) on  qubit %zu is %zu\n", i, flagCirc->len_qcirc[i]);
		printf("-------------\n");

        for(size_t j = 0; j < flagCirc->len_qcirc[i]; j++){

            printf("%s %zu ", flagCirc->ins_id[flagCirc->qcirc[i][j]].name, flagCirc->ins_id[flagCirc->qcirc[i][j]].op[0]);

            if(strcmp(flagCirc->ins_id[flagCirc->qcirc[i][j]].name, "cnot") == 0 
               || strcmp(flagCirc->ins_id[flagCirc->qcirc[i][j]].name, "ch") == 0){

                printf("%zu\n", flagCirc->ins_id[flagCirc->qcirc[i][j]].op[1]);

            }else{

                printf("\n");

            }

        }
	}

    printf("Ancilla is qubit %zu\n", flagCirc->ancillaInd);
    printf("-------------\n");

    return;
}
