#include <stdio.h>
#include <stdlib.h>
#include <regex.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <time.h>  

#include "omp.h"
#include "parser.h"
#include "error_propagation.h"

/* Function declaration */
size_t identify_ancilla(flag_circ flagCirc);
void insert_flag(size_t first_pos, size_t second_pos, flag_circ * flagCirc);
size_t calculate_syndrome(error_t error_result[], flag_circ * flagCirc);
void build_error_dict(flag_circ * flagCirc, syndrome_ind * synInd);
void update_error_dict(size_t i, error_t ET, flag_circ * flagCirc, syndrome_ind * synInd);
void print_error_dict_info(size_t maxE, syndrome_ind synInd, flag_circ flagCirc);
bool test_result(size_t maxE, flag_circ flagCirc, error_t error_result[]);
bool random_flag_test(flag_circ * flagCirc, size_t maxE);
bool test_flag_circuit(flag_circ newCirc, size_t maxE);
bool simple_loop_check(size_t inds[16], size_t num_data_ch);
size_t random_num(size_t ancLen);

/* Main function */
int main(int argc, char *argv[]){

    /* Set seed for random number generation */
    srand(time(NULL));    

    /* Define variables */
    size_t maxError=INVALID_SIZE;
    char fileName[MAX_FILENAME_SIZE];

    flag_circ flagCirc={.qcirc={{INVALID_SIZE}},
                        .len_qcirc={0},
                        .dataQubitNum=0,
                        .flagQubitNum=0,
                        .currentQubitNum=0,
                        .ancillaInd=INVALID_SIZE,
                        .len_ins=0
                       };

    raw_error rawError={.len={0}};
    error_t error_result[MAX_QUBIT_NUM];

    /* Parse input */
    handle_input(argc, argv, &maxError, &flagCirc, fileName);
    parse_qasm(fileName, &flagCirc);
    flagCirc.currentQubitNum = flagCirc.dataQubitNum + 1;
    flagCirc.ancillaInd = identify_ancilla(flagCirc); // Find the qubit prepared in |+> and measured in X basis. Return the first match.

    bool res = false;
    while(true){
       random_flag_test(&flagCirc, maxError);
    }

    return 0;
}

/* Function implementation */
size_t identify_ancilla(flag_circ flagCirc){
    
    for(size_t i = 0; i < flagCirc.dataQubitNum + 1; i++){
        if(flagCirc.len_qcirc[i] == 0) continue;
        if(strcmp(flagCirc.ins_id[flagCirc.qcirc[i][0]].name, "prepx") == 0 &&
           strcmp(flagCirc.ins_id[flagCirc.qcirc[i][flagCirc.len_qcirc[i]-1]].name, "measx") == 0){
           return i; 
        }
    }
    
    printf("Didn't find an ancilla, please make sure this is a logical Hadamard measurement circuit\n");
    exit(0);
}

void insert_flag(size_t first_pos, size_t second_pos, flag_circ * flagCirc){
    
    if(first_pos >= second_pos) return;

    size_t ctrl = flagCirc->ancillaInd, targ = flagCirc->currentQubitNum;
    size_t prepz_id = flagCirc->len_ins, cnot1_id = flagCirc->len_ins + 1, cnot2_id = flagCirc->len_ins + 2, measz_id = flagCirc->len_ins + 3;

    flagCirc->len_ins += 4;
    flagCirc->currentQubitNum += 1;

    flagCirc->len_qcirc[ctrl] += 2;

    flagCirc->qcirc[targ][flagCirc->len_qcirc[targ]] = prepz_id;
    flagCirc->len_qcirc[targ] += 1;
    flagCirc->qcirc[targ][flagCirc->len_qcirc[targ]] = cnot1_id;
    flagCirc->len_qcirc[targ] += 1;
    flagCirc->qcirc[targ][flagCirc->len_qcirc[targ]] = cnot2_id;
    flagCirc->len_qcirc[targ] += 1;
    flagCirc->qcirc[targ][flagCirc->len_qcirc[targ]] = measz_id;
    flagCirc->len_qcirc[targ] += 1;

    strcpy(flagCirc->ins_id[prepz_id].name, "prepz");
    strcpy(flagCirc->ins_id[cnot1_id].name, "cnot");
    strcpy(flagCirc->ins_id[cnot2_id].name, "cnot");
    strcpy(flagCirc->ins_id[measz_id].name, "measz");

    flagCirc->ins_id[prepz_id].op[0] = targ;
    flagCirc->ins_id[prepz_id].qcirc_ind[0] = 0;

    flagCirc->ins_id[cnot1_id].op[0] = ctrl;
    flagCirc->ins_id[cnot1_id].op[1] = targ;
    flagCirc->ins_id[cnot1_id].qcirc_ind[0] = first_pos;
    flagCirc->ins_id[cnot1_id].qcirc_ind[1] = 1;

    flagCirc->ins_id[cnot2_id].op[0] = ctrl;
    flagCirc->ins_id[cnot2_id].op[1] = targ;
    flagCirc->ins_id[cnot2_id].qcirc_ind[0] = second_pos + 1;
    flagCirc->ins_id[cnot2_id].qcirc_ind[1] = 2;

    flagCirc->ins_id[measz_id].op[0] = targ;
    flagCirc->ins_id[measz_id].qcirc_ind[0] = 3;
    
    for(size_t i = flagCirc->len_qcirc[ctrl] -3; i >= second_pos; i--){

       size_t gate_id = flagCirc->qcirc[ctrl][i];

       flagCirc->qcirc[ctrl][i+2] = gate_id;

       if(flagCirc->ins_id[gate_id].op[0] == ctrl)
           flagCirc->ins_id[gate_id].qcirc_ind[0] += 2;
       if(flagCirc->ins_id[gate_id].op[1] == ctrl)
           flagCirc->ins_id[gate_id].qcirc_ind[1] += 2;
    }
    
    flagCirc->qcirc[ctrl][second_pos+1] = cnot2_id;

    for(size_t i = second_pos - 1; i >= first_pos; i--){

       size_t gate_id = flagCirc->qcirc[ctrl][i];

       flagCirc->qcirc[ctrl][i+1] = gate_id;
       if(flagCirc->ins_id[gate_id].op[0] == ctrl){
           flagCirc->ins_id[gate_id].qcirc_ind[0] = flagCirc->ins_id[gate_id].qcirc_ind[0] + 1;
       }
       else if(flagCirc->ins_id[gate_id].op[1] == ctrl){
           flagCirc->ins_id[gate_id].qcirc_ind[1] = flagCirc->ins_id[gate_id].qcirc_ind[1] + 1;
       }
         
    }

    flagCirc->qcirc[ctrl][first_pos] = cnot1_id;
}

/* This function assumes only cnot, ch, in our case, preparation and measurement error are equivalent to cnot error */
void build_error_dict(flag_circ * flagCirc, syndrome_ind * synInd){

    for(size_t i = 0; i < flagCirc->len_ins; i++){ 

        if(strcmp(flagCirc->ins_id[i].name, "cnot")==0 || strcmp(flagCirc->ins_id[i].name, "ch")==0){
            update_error_dict(i, IX, flagCirc, synInd);
            update_error_dict(i, XI, flagCirc, synInd);
            update_error_dict(i, XX, flagCirc, synInd);
        }
    }
    return;
}

void update_error_dict(size_t i, error_t ET, flag_circ * flagCirc, syndrome_ind * synInd){

    error_t error_result[MAX_QUBIT_NUM]={I};
    raw_error rawError={.len={0}};

    error_gen(i, ET, flagCirc, &rawError); 
    raw_to_final_error_result(flagCirc->currentQubitNum, &rawError, error_result);

    size_t ext_flag = 1;

    for(size_t j = 0; j < synInd->len; j++){
        size_t eq_flag = 1;
        for(size_t k = 0; k < flagCirc->currentQubitNum; k++){
            
            if(error_result[k]!=synInd->synHT[j].syndrome[k]){
                eq_flag = -1;
                break;
            }
        }
        if(eq_flag == 1){
            ext_flag = -1; 
            synInd->synHT[j].loc[synInd->synHT[j].loc_len] = i;
            synInd->synHT[j].error_type[synInd->synHT[j].loc_len] = ET;
            synInd->synHT[j].loc_len += 1;
            break;
        }
    }

    if(ext_flag == 1){
        for(size_t k = 0; k < flagCirc->currentQubitNum; k++){
            synInd->synHT[synInd->len].syndrome[k] = error_result[k];
        }
        synInd->synHT[synInd->len].loc[synInd->synHT[synInd->len].loc_len] = i;
        synInd->synHT[synInd->len].error_type[synInd->synHT[synInd->len].loc_len] = ET;
        synInd->synHT[synInd->len].loc_len += 1;
        synInd->len += 1;
    }
    return;
}

void print_error_dict_info(size_t maxE, syndrome_ind synInd, flag_circ flagCirc){

    printf("There are in total %zu syndromes.\n", synInd.len);
    printf("----------------------");
    
    for(size_t i = 0 ; i< synInd.len;i++){
        printf("\nSyndrome %zu is:\n", i);
        for(size_t j =0; j < flagCirc.currentQubitNum; j++){
            char error_string[10];
            error_to_string(synInd.synHT[i].syndrome[j], error_string);
            printf("%s  ", error_string);
        }
        printf("There are these errors that produce this syndrome.\n");
        for(size_t j =0; j < synInd.synHT[i].loc_len; j++){
            char error_string[3];
            error_to_string(synInd.synHT[i].error_type[j], error_string);
            printf("gate %s %zu %zu %s error;", 
                   flagCirc.ins_id[synInd.synHT[i].loc[j]].name, 
                   flagCirc.ins_id[synInd.synHT[i].loc[j]].op[0],
                   flagCirc.ins_id[synInd.synHT[i].loc[j]].op[1],
                   error_string);
        }
        
        bool res;
        res = test_result(maxE, flagCirc, synInd.synHT[i].syndrome);
        printf("This syndrome is good? %s", res ? "true" : "false");
    }
}

bool test_result(size_t maxE, flag_circ flagCirc, error_t error_result[]){

    bool flagged = false;
    for(size_t i = flagCirc.dataQubitNum + 1; i < flagCirc.currentQubitNum; i++){
        if(error_result[i] != I){
            flagged = true; 
        }  
    }
    
    size_t error_support = 0, flipped_error_support = 0;
    for(size_t i = 0; i < flagCirc.dataQubitNum; i++){
        if(error_result[i] != I){
            error_support += 1; 
        } 
        if(combine_error(H, error_result[i])!=I){
            flipped_error_support += 1; 

        }
    }

    if((error_support > maxE) && (flipped_error_support > maxE) && !flagged){
        return false;        
    }

    return true;
}

bool random_flag_test(flag_circ * flagCirc, size_t maxE){
    
    flag_circ newCirc;
    memcpy((void *) &newCirc, (void *) flagCirc, sizeof(newCirc));

    for(size_t i = 0; i < newCirc.flagQubitNum; i++){

        size_t ancInd = newCirc.ancillaInd;
        size_t ancLen = newCirc.len_qcirc[ancInd];
        
        size_t first = 0, second = 0;
        while(first == second){
            first = random_num(ancLen);
            second = random_num(ancLen);
        }

        if(first > second){
            size_t temp;
            temp = first;
            first = second;
            second = temp;
        }

        insert_flag(first, second, &newCirc);
    }
    
    bool res;
    res = test_flag_circuit(newCirc, maxE);
    if(res) print_summary(maxE, &newCirc);

    return res;
}

bool test_flag_circuit(flag_circ newCirc, size_t maxE){

    /* Syndrome hash table */
    syndrome_ind synInd={.len=0};

    /* Initialization of the hash table */
    for(size_t j = 0; j < MAX_INSTRUCTION_NUM * 3; j++){
        synInd.synHT[j].loc_len = 0;
    }
    
    build_error_dict(&newCirc, &synInd);

    // print_summary(maxE, &newCirc);

    // print_error_dict_info(maxE, synInd, newCirc);
    
    if(maxE == 2){
        for(size_t i = 0; i < synInd.len - 1; i++){
            for(size_t j = i + 1; j < synInd.len; j++){
                error_t error_result[newCirc.currentQubitNum];

                for(size_t m = 0;  m < newCirc.currentQubitNum; m++){

                    error_result[m] = synInd.synHT[i].syndrome[m];
                    error_result[m] = combine_error(synInd.synHT[j].syndrome[m], 
                                                    error_result[m]);
                }

                bool res;
                res = test_result(maxE, newCirc, error_result);
                if(!res) {
                    return false; 
                }
            }
        } 
    }

    else if(maxE == 3){
        for(size_t i = 0; i < synInd.len - 2; i++){
            for(size_t j = i + 1; j < synInd.len - 1; j++){
                for(size_t k = j + 1; k < synInd.len; k++){
                    error_t error_result[newCirc.currentQubitNum];
                    for(size_t m = 0;  m < newCirc.currentQubitNum; m++){

                        error_result[m] = synInd.synHT[i].syndrome[m];
                        error_result[m] = combine_error(synInd.synHT[j].syndrome[m], 
                                                        error_result[m]);
                        error_result[m] = combine_error(synInd.synHT[k].syndrome[m], 
                                                        error_result[m]);
                        
                    }
                    bool res;
                    res = test_result(maxE, newCirc, error_result);
                    if(!res) {return false;}
                }
            }
        }
    }

    else if(maxE == 1){
        for(size_t i = 0; i < synInd.len; i++){
            error_t error_result[newCirc.currentQubitNum];
            for(size_t m = 0;  m < newCirc.currentQubitNum; m++){
                error_result[m] = synInd.synHT[i].syndrome[m];
            }
            bool res;
            res = test_result(maxE, newCirc, error_result);
            if(!res) {return false;}
        }
    
    }
    else{
        printf("Invalid max error number.\n");
        exit(0);
    }
    
    return true;
}

bool simple_loop_check(size_t inds[16], size_t num_data_ch){
   /* i11 i12 i21 i22 i31 i32 i41 i42 ...*/

    if(inds[1] <= inds[0]) return false;
    if(inds[3] <= inds[2]) return false;
    if(inds[5] <= inds[4]) return false;
    if(inds[7] <= inds[6]) return false;
    if(inds[9] <= inds[8]) return false;
    if(inds[11] <= inds[10]) return false;
    if(inds[13] <= inds[12]) return false;
    if(inds[15] <= inds[14]) return false;

    if(inds[2] < inds[0]) return false;
    if(inds[4] < inds[2]) return false;
    if(inds[6] < inds[4]) return false;
    if(inds[8] < inds[6]) return false;
    if(inds[10] < inds[8]) return false;
    if(inds[12] < inds[10]) return false;
    if(inds[14] < inds[12]) return false;
    
    size_t max = inds[0];

    for(size_t i =1; i < 16; i+=2){
        if(max < inds[i])
            max = inds[i];
    }

    if(max < num_data_ch + 14 - 3) return false;

    return true;
}

size_t random_num(size_t ancLen){
    return (size_t) (rand() % (ancLen - 1)) + 1;
}
