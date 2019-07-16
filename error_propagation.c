#include <stdio.h>
#include <stdlib.h>
#include <regex.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include "error_propagation.h"

void clifford_propagation(error_t error_pauli, size_t qbit_index, size_t depth,
                          flag_circ * flagCirc, raw_error * rawError){
	
	size_t pos = flagCirc->qcirc[qbit_index][depth];
	size_t ctrl, targ, ctrl_ind, targ_ind;
    
    /* Propgate here through the last gate. Unlike measz, measx, keep all the errors*/
    if(depth == flagCirc->len_qcirc[qbit_index]){

        rawError->result[qbit_index][rawError->len[qbit_index]] = error_pauli;
        rawError->len[qbit_index] += 1;
		return;

    }
    
    /* Measx */
    if(strcmp(flagCirc->ins_id[pos].name, "measx") == 0){

		if(error_pauli == Y || error_pauli == Z){
			rawError->result[qbit_index][rawError->len[qbit_index]] = Z;
			rawError->len[qbit_index] += 1;
		}

		else if(error_pauli == A){
			rawError->result[qbit_index][rawError->len[qbit_index]] = A;
			rawError->len[qbit_index] += 1;
		}

		else if(error_pauli == H || error_pauli == HX){
			rawError->result[qbit_index][rawError->len[qbit_index]] = H;
			rawError->len[qbit_index] += 1;
		}

        else if (error_pauli == HZ || error_pauli == HY){
			rawError->result[qbit_index][rawError->len[qbit_index]] = HZ;
			rawError->len[qbit_index] += 1;
        }

		return;
    }

    /* Measz*/
    if(strcmp(flagCirc->ins_id[pos].name, "measz") == 0){

		if(error_pauli == Y | error_pauli == X){
			rawError->result[qbit_index][rawError->len[qbit_index]] = X;
			rawError->len[qbit_index] += 1;
		}
		else if(error_pauli == A){
			rawError->result[qbit_index][rawError->len[qbit_index]] = A;
			rawError->len[qbit_index] += 1;
		}
		else if(error_pauli == H || error_pauli == HZ){
			rawError->result[qbit_index][rawError->len[qbit_index]] = H;
			rawError->len[qbit_index] += 1;
		}

        else if (error_pauli == HY || error_pauli == HX){
			rawError->result[qbit_index][rawError->len[qbit_index]] = HX;
			rawError->len[qbit_index] += 1;
        }

		return;
    }
    
    /* For CNOT */
	if((pos > 0 || pos == 0) && strcmp(flagCirc->ins_id[pos].name, "cnot") == 0){

        size_t ctrl_depth, targ_depth;

		ctrl = flagCirc->ins_id[pos].op[0];
		targ = flagCirc->ins_id[pos].op[1];
		
		ctrl_ind = flagCirc->ins_id[pos].qcirc_ind[0];
		targ_ind = flagCirc->ins_id[pos].qcirc_ind[1];

        ctrl_depth = flagCirc->len_qcirc[ctrl] - 1;
        targ_depth = flagCirc->len_qcirc[targ] - 1;

        /* If an X error pass through a control bit or a Z error pass through a target bit,
         * the error has to propagate to the other bit */
		if((ctrl == qbit_index && error_pauli == X) || (targ == qbit_index && error_pauli == Z)){
            clifford_propagation(error_pauli, ctrl, ctrl_ind + 1, flagCirc, rawError);
			clifford_propagation(error_pauli, targ, targ_ind + 1, flagCirc, rawError);
			return;
		}

        /* If a Z error pass through a control bit, the error commutes */
		if(ctrl == qbit_index && error_pauli == Z){
			clifford_propagation(error_pauli, ctrl, ctrl_ind + 1, flagCirc, rawError);
			return;
		}
		
        /* If an X error pass through a target bit, the error commutes */
		if(targ == qbit_index && error_pauli == X){
			clifford_propagation(error_pauli, targ, targ_ind + 1, flagCirc, rawError);
			return;
		}
        
        /* Y error: YI --> YX; IY --> ZY */
		if(error_pauli == Y){
			if(ctrl == qbit_index){
				clifford_propagation(Y, ctrl, ctrl_ind + 1, flagCirc, rawError);
				clifford_propagation(X, targ, targ_ind + 1, flagCirc, rawError);
				return;
			}
			if(targ == qbit_index){
				clifford_propagation(Z, ctrl, ctrl_ind + 1, flagCirc, rawError);
				clifford_propagation(Y, targ, targ_ind + 1, flagCirc, rawError);
				return;
			}
		}

        if(error_pauli == H || error_pauli == A || error_pauli == HX || error_pauli == HY || error_pauli == HZ){
            clifford_propagation(A, ctrl, ctrl_ind + 1, flagCirc, rawError);
            clifford_propagation(A, targ, targ_ind + 1, flagCirc, rawError);
            return;
        }

	}	
    
    /* For controlled-Hadamard */
	if((pos > 0 || pos == 0 ) && strcmp(flagCirc->ins_id[pos].name, "ch") == 0){

		ctrl = flagCirc->ins_id[pos].op[0];
		targ = flagCirc->ins_id[pos].op[1];
		
		ctrl_ind = flagCirc->ins_id[pos].qcirc_ind[0];
		targ_ind = flagCirc->ins_id[pos].qcirc_ind[1];
    
        /* XI --> XH */
		if(ctrl == qbit_index && error_pauli == X){
			clifford_propagation(X, ctrl, ctrl_ind + 1, flagCirc, rawError);
			clifford_propagation(H, targ, targ_ind + 1, flagCirc, rawError);
			return;
		}
        
        /* IX --> AA
         * Reason: (P0 (X) I + P1 (X) H) * (I (X) X) = P0 (X) X + P1 (X) HZ 
         * * which is half X and half Z, we denote it as A(bnormal).
         * */
		if(targ == qbit_index && error_pauli == X){
			clifford_propagation(A, ctrl, ctrl_ind + 1, flagCirc, rawError);
			clifford_propagation(A, targ, targ_ind + 1, flagCirc, rawError);
			return;
		}

        /* ZI --> ZI  */
		if(ctrl == qbit_index && error_pauli == Z){
			clifford_propagation(error_pauli, ctrl, ctrl_ind + 1, flagCirc, rawError);
			return;
		}

        /* IZ --> AA
         * Reason: 
         * which is half X and half Z, we denote it as A(bnormal).
         * */
		if(targ == qbit_index && error_pauli == Z){
			clifford_propagation(A, targ, targ_ind + 1, flagCirc, rawError);
			clifford_propagation(A, ctrl, ctrl_ind + 1, flagCirc, rawError);
			return;
		}

        /* Y errors */
		if(error_pauli == Y){
			if(ctrl == qbit_index){
				clifford_propagation(Y, ctrl, ctrl_ind + 1, flagCirc, rawError);
				clifford_propagation(A, targ, targ_ind + 1, flagCirc, rawError);
				return;
			}
			if(targ == qbit_index){
				clifford_propagation(Z, ctrl, ctrl_ind + 1, flagCirc, rawError);
				clifford_propagation(Y, targ, targ_ind + 1, flagCirc, rawError);
				return;
			}
		}

        if(error_pauli == A){
            clifford_propagation(A, ctrl, ctrl_ind + 1, flagCirc, rawError);
            clifford_propagation(A, targ, targ_ind + 1, flagCirc, rawError);
            return;
        }

        if(error_pauli == H){
			if(ctrl == qbit_index){
				clifford_propagation(A, ctrl, ctrl_ind + 1, flagCirc, rawError);
				clifford_propagation(A, targ, targ_ind + 1, flagCirc, rawError);
				return;
			}
			if(targ == qbit_index){
				clifford_propagation(H, targ, targ_ind + 1, flagCirc, rawError);
				return;
			}
        }

        if(error_pauli == HX){
            clifford_propagation(A, ctrl, ctrl_ind + 1, flagCirc, rawError);
            clifford_propagation(A, targ, targ_ind + 1, flagCirc, rawError);
            return;
        }

        if(error_pauli == HY){
			if(ctrl == qbit_index){
				clifford_propagation(A, ctrl, ctrl_ind + 1, flagCirc, rawError);
				clifford_propagation(A, targ, targ_ind + 1, flagCirc, rawError);
				return;
			}
			if(targ == qbit_index){
				clifford_propagation(Z, ctrl, ctrl_ind + 1, flagCirc, rawError);
				clifford_propagation(HY, targ, targ_ind + 1, flagCirc, rawError);
				return;
			}
        }

        if(error_pauli == HZ){
			if(ctrl == qbit_index){
				clifford_propagation(A, ctrl, ctrl_ind + 1, flagCirc, rawError);
				clifford_propagation(A, targ, targ_ind + 1, flagCirc, rawError);
				return;
			}
			if(targ == qbit_index){
				clifford_propagation(Z, ctrl, ctrl_ind + 1, flagCirc, rawError);
				clifford_propagation(HY, targ, targ_ind + 1, flagCirc, rawError);
				return;
			}
        }

	}
    
    /* Hadamard */
	if((pos > 0 || pos == 0 ) && strcmp(flagCirc->ins_id[pos].name, "h") == 0){

        switch(error_pauli){
        
            case X: 
                clifford_propagation(Z, qbit_index, depth + 1, flagCirc, rawError);
                break;

            case Y:
                clifford_propagation(Y, qbit_index, depth + 1, flagCirc, rawError);
                break;
            
            case Z:
                clifford_propagation(X, qbit_index, depth + 1, flagCirc, rawError);
                break;

            case H:
                clifford_propagation(H, qbit_index, depth + 1, flagCirc, rawError);
                break;

            case A:
                clifford_propagation(A, qbit_index, depth + 1, flagCirc, rawError);
                break;

            case HX:
                clifford_propagation(HZ, qbit_index, depth + 1, flagCirc, rawError);
                break;

            case HY:
                clifford_propagation(HY, qbit_index, depth + 1, flagCirc, rawError);
                break;

            case HZ:
                clifford_propagation(HX, qbit_index, depth + 1, flagCirc, rawError);

        }

        return;
    }

    /* S gate */
	if((pos > 0 || pos == 0 ) && strcmp(flagCirc->ins_id[pos].name, "s") == 0){
        switch(error_pauli){
        
            case X: 
                clifford_propagation(Y, qbit_index, depth + 1, flagCirc, rawError);
                break;

            case Y:
                clifford_propagation(X, qbit_index, depth + 1, flagCirc, rawError);
                break;
            
            case Z:
                clifford_propagation(Z, qbit_index, depth + 1, flagCirc, rawError);
                break;

            case H:
                clifford_propagation(A, qbit_index, depth + 1, flagCirc, rawError);
                break;

            case A:
                clifford_propagation(A, qbit_index, depth + 1, flagCirc, rawError);
                break;

            case HX:
                clifford_propagation(A, qbit_index, depth + 1, flagCirc, rawError);
                break;

            case HY:
                clifford_propagation(A, qbit_index, depth + 1, flagCirc, rawError);
                break;

            case HZ:
                clifford_propagation(A, qbit_index, depth + 1, flagCirc, rawError);

        }
        return;
    } 

    /* Paulis */
	if((pos > 0 || pos == 0 ) && ((strcmp(flagCirc->ins_id[pos].name, "x") == 0 || strcmp(flagCirc->ins_id[pos].name, "z") == 0)
       || strcmp(flagCirc->ins_id[pos].name, "y") == 0)){
        switch(error_pauli){
        
            case X: 
                clifford_propagation(X, qbit_index, depth + 1, flagCirc, rawError);
                break;

            case Y:
                clifford_propagation(Y, qbit_index, depth + 1, flagCirc, rawError);
                break;
            
            case Z:
                clifford_propagation(Z, qbit_index, depth + 1, flagCirc, rawError);
                break;

            case A:
                clifford_propagation(A, qbit_index, depth + 1, flagCirc, rawError);
        }
	}

	if((pos > 0 || pos == 0 ) && (strcmp(flagCirc->ins_id[pos].name, "x") == 0)){
        switch(error_pauli){
            case H:
                clifford_propagation(HY, qbit_index, depth + 1, flagCirc, rawError);
                break;
            case HX:
                clifford_propagation(HZ, qbit_index, depth + 1, flagCirc, rawError);
                break;
            case HY:
                clifford_propagation(H, qbit_index, depth + 1, flagCirc, rawError);
                break;
            case HZ:
                clifford_propagation(HX, qbit_index, depth + 1, flagCirc, rawError);
        }
        return;
    }

	if((pos > 0 || pos == 0 ) && (strcmp(flagCirc->ins_id[pos].name, "y") == 0)){
        switch(error_pauli){
            case H:
                clifford_propagation(H, qbit_index, depth + 1, flagCirc, rawError);
                break;
            case HX:
                clifford_propagation(HX, qbit_index, depth + 1, flagCirc, rawError);
                break;
            case HY:
                clifford_propagation(HY, qbit_index, depth + 1, flagCirc, rawError);
                break;
            case HZ:
                clifford_propagation(HZ, qbit_index, depth + 1, flagCirc, rawError);
        }
        return;
    }

	if((pos > 0 || pos == 0 ) && (strcmp(flagCirc->ins_id[pos].name, "z") == 0)){
        switch(error_pauli){
            case H:
                clifford_propagation(HY, qbit_index, depth + 1, flagCirc, rawError);
                break;
            case HX:
                clifford_propagation(HZ, qbit_index, depth + 1, flagCirc, rawError);
                break;
            case HY:
                clifford_propagation(HX, qbit_index, depth + 1, flagCirc, rawError);
                break;
            case HZ:
                clifford_propagation(H, qbit_index, depth + 1, flagCirc, rawError);
        }
        return;
    }

	return;
}

/* Error generation */
void error_gen(size_t gate_id, error_t error_type, flag_circ * flagCirc, raw_error * rawError){

	if(strcmp(flagCirc->ins_id[gate_id].name, "cnot") == 0 ||strcmp(flagCirc->ins_id[gate_id].name, "ch") == 0){
    
        size_t ctrl, targ, ctrl_ind, targ_ind;

		ctrl = flagCirc->ins_id[gate_id].op[0];
		targ = flagCirc->ins_id[gate_id].op[1];
		
		ctrl_ind = flagCirc->ins_id[gate_id].qcirc_ind[0];
		targ_ind = flagCirc->ins_id[gate_id].qcirc_ind[1];
        
        char error_string[15];

        switch(error_type){
            case IX:
                clifford_propagation(X, targ, targ_ind + 1, flagCirc, rawError);
                break;
            case XI:
                clifford_propagation(X, ctrl, ctrl_ind + 1, flagCirc, rawError);
                break;
            case XX:
                clifford_propagation(X, ctrl, ctrl_ind + 1, flagCirc, rawError);
                clifford_propagation(X, targ, targ_ind + 1, flagCirc, rawError);
                break;
            case XY:
                clifford_propagation(X, ctrl, ctrl_ind + 1, flagCirc, rawError);
                clifford_propagation(Y, targ, targ_ind + 1, flagCirc, rawError);
                break;
            case XZ:
                clifford_propagation(X, ctrl, ctrl_ind + 1, flagCirc, rawError);
                clifford_propagation(Z, targ, targ_ind + 1, flagCirc, rawError);
                break;
            case IZ:
                clifford_propagation(Z, targ, targ_ind + 1, flagCirc, rawError);
                break;
            case ZI:
                clifford_propagation(Z, ctrl, ctrl_ind + 1, flagCirc, rawError);
                break;
            case ZZ:
                clifford_propagation(Z, ctrl, ctrl_ind + 1, flagCirc, rawError);
                clifford_propagation(Z, targ, targ_ind + 1, flagCirc, rawError);
                break;
            case ZX:
                clifford_propagation(Z, ctrl, ctrl_ind + 1, flagCirc, rawError);
                clifford_propagation(X, targ, targ_ind + 1, flagCirc, rawError);
                break;
            case ZY:
                clifford_propagation(Z, ctrl, ctrl_ind + 1, flagCirc, rawError);
                clifford_propagation(Y, targ, targ_ind + 1, flagCirc, rawError);
                break;
            case IY:
                clifford_propagation(Y, targ, targ_ind + 1, flagCirc, rawError);
                break;
            case YI:
                clifford_propagation(Y, ctrl, ctrl_ind + 1, flagCirc, rawError);
                break;
            case YY:
                clifford_propagation(Y, ctrl, ctrl_ind + 1, flagCirc, rawError);
                clifford_propagation(Y, targ, targ_ind + 1, flagCirc, rawError);
            case YX:
                clifford_propagation(Y, ctrl, ctrl_ind + 1, flagCirc, rawError);
                clifford_propagation(X, targ, targ_ind + 1, flagCirc, rawError);
                break;
            case YZ:
                clifford_propagation(Y, ctrl, ctrl_ind + 1, flagCirc, rawError);
                clifford_propagation(Z, targ, targ_ind + 1, flagCirc, rawError);
                break;

            default:
                error_to_string(error_type, error_string);
                printf("Invalid error type %s on %s on qubit %zu(%zu th gate, including prep) and %zu(%zu th gate, including prep). We only have Pauli error generation. \n", error_string, flagCirc->ins_id[gate_id].name, ctrl, ctrl_ind + 1, targ, targ_ind + 1);
                exit(0);
        }
        return;

	}

	if(strcmp(flagCirc->ins_id[gate_id].name, "x") == 0 || strcmp(flagCirc->ins_id[gate_id].name, "y") == 0 ||
       strcmp(flagCirc->ins_id[gate_id].name, "z") == 0 || strcmp(flagCirc->ins_id[gate_id].name, "h") == 0 ||
       strcmp(flagCirc->ins_id[gate_id].name, "s") == 0){

        size_t qbit_index, depth;

        qbit_index = flagCirc->ins_id[gate_id].op[0];
        depth = flagCirc->ins_id[gate_id].qcirc_ind[0];

        if(depth == (flagCirc->len_qcirc[qbit_index] - 1)){

            if(error_type == Y | error_type == X){
                rawError->result[qbit_index][rawError->len[qbit_index]] = X;
                rawError->len[qbit_index] += 1;
            }
            return;
        }
        
        char error_string[15];
        switch(error_type){ 
            case X:
				clifford_propagation(X, qbit_index, depth + 1, flagCirc, rawError);
				break;
			case Y:
				clifford_propagation(Y, qbit_index, depth + 1, flagCirc, rawError);
				break;
			case Z:
				clifford_propagation(Z, qbit_index, depth + 1, flagCirc, rawError);
                break;
            default:
                error_to_string(error_type, error_string);
                printf("Invalid error type %s on the %zu th gate on qubit %zu. We only have Pauli error channel generation.\n", error_string, depth + 1, qbit_index);
                exit(0);
        }
        return;
    }

	if(strcmp(flagCirc->ins_id[gate_id].name, "prepx") == 0){

        size_t qbit_index, depth;

        qbit_index = flagCirc->ins_id[gate_id].op[0];
        depth = flagCirc->ins_id[gate_id].qcirc_ind[0];

        if(depth == (flagCirc->len_qcirc[qbit_index] - 1)){

            if(error_type == Y | error_type == Z){
                rawError->result[qbit_index][rawError->len[qbit_index]] = Z;
                rawError->len[qbit_index] += 1;
            }
            return;
        }

        char error_string[15];    
        switch(error_type){ 
			case Y:
				clifford_propagation(Z, qbit_index, depth + 1, flagCirc, rawError);
				break;
			case Z:
				clifford_propagation(Z, qbit_index, depth + 1, flagCirc, rawError);
                break;
            default:
                error_to_string(error_type, error_string);
                printf("Invalid error type %s on %zu th gate(including prep) on qubit %zu. After prepx, we only have error Y, Z, they result in Z error.\n", error_string, depth + 1, qbit_index);
                exit(0);
        }
        return;
    
    }

	if(strcmp(flagCirc->ins_id[gate_id].name, "prepz") == 0){

        size_t qbit_index, depth;

        qbit_index = flagCirc->ins_id[gate_id].op[0];
        depth = flagCirc->ins_id[gate_id].qcirc_ind[0];

        if(depth == (flagCirc->len_qcirc[qbit_index] - 1)){

            if(error_type == X | error_type == Y){
                rawError->result[qbit_index][rawError->len[qbit_index]] = X;
                rawError->len[qbit_index] += 1;
            }
            return;
        }
        
        char error_string[15];
        switch(error_type){ 
			case X:
				clifford_propagation(X, qbit_index, depth + 1, flagCirc, rawError);
				break;
			case Y:
				clifford_propagation(X, qbit_index, depth + 1, flagCirc, rawError);
                break;
            default:
                printf("Invalid error type %s on %zu th gate(including prep) on qubit %zu. After prepz, we only have error Y, X, they result in X error.\n", error_string, depth + 1, qbit_index);
                exit(0);
        }
        return;
    
    }
}

void print_raw_error_result(size_t qubitNum, raw_error * rawError){

    printf("The errors on each qubit is as follow:\n");
    printf("--------------------------------------");

    for (size_t i = 0; i < qubitNum; i++){
        printf("\nOn qubit %zu:", i);
        for(size_t j = 0; j < rawError->len[i]; j++){

            char error_string[15];
            error_to_string(rawError->result[i][j], error_string);
            printf("%s ", error_string);
        }
    }
}

void error_to_string(error_t raw, char error_string[]){
    switch(raw){
        case X:
            strcpy(error_string, "X\0");
            break;
        case Y:
            strcpy(error_string, "Y\0");
            break;
        case Z:
            strcpy(error_string, "Z\0");
            break;
        case I:
            strcpy(error_string, "I\0");
            break;
        case H:
            strcpy(error_string, "H\0");
            break;
        case A:
            strcpy(error_string, "A\0");
            break;
        case HX:
            strcpy(error_string, "HX\0");
            break;
        case HY:
            strcpy(error_string, "HY\0");
            break;
        case HZ:
            strcpy(error_string, "HZ\0");
            break;
        case IX:
            strcpy(error_string, "IX\0");
            break;
        case IY:
            strcpy(error_string, "IY\0");
            break;
        case IZ:
            strcpy(error_string, "IZ\0");
            break;
        case XX:
            strcpy(error_string, "XX\0");
            break;
        case XY:
            strcpy(error_string, "XY\0");
            break;
        case XZ:
            strcpy(error_string, "XZ\0");
            break;
        case XI:
            strcpy(error_string, "XI\0");
            break;
        case YX:
            strcpy(error_string, "YX\0");
            break;
        case YY:
            strcpy(error_string, "YY\0");
            break;
        case YZ:
            strcpy(error_string, "YZ\0");
            break;
        case YI:
            strcpy(error_string, "YI\0");
            break;
        case ZX:
            strcpy(error_string, "ZX\0");
            break;
        case ZY:
            strcpy(error_string, "ZY\0");
            break;
        case ZZ:
            strcpy(error_string, "ZZ\0");
            break;
        case ZI:
            strcpy(error_string, "ZI\0");
            break;
        case N:
            strcpy(error_string, "Not specified\0");
    } 
    return;
}

void raw_to_final_error_result(size_t qubitNum, raw_error * rawError, error_t error_result[]){
    if(qubitNum < 1){
        printf("Measuring no qubits. Error convertin raw error result to final error result.\n"); 
    }
    for(size_t i  = 0; i < qubitNum; i++){
        if(rawError->len[i] == 0) error_result[i] = N;
        
        error_t rt_val = rawError->result[i][0];
        for(size_t j = 1; j < rawError->len[i]; j++){
            rt_val = combine_error(rawError->result[i][j], rt_val);
        }
        error_result[i] = rt_val;
    }

    return;
}

error_t combine_error(error_t first, error_t second){
    if(first == I) return second;
    if(second == I) return first;
    if(first == A || second == A) return A;
    if(first == N || second == N) return N;

    if(first == X){
        if(second == X) return I;
        if(second == Y) return Z;
        if(second == Z) return Y;
        if(second == H) return HZ;
        if(second == HX) return HY;
        if(second == HZ) return H;
        if(second == HY) return HX;
    }

    if(first == Y){
    
        if(second == X) return Z;
        if(second == Y) return I;
        if(second == Z) return X;
        if(second == H) return HY;
        if(second == HX) return HZ;
        if(second == HZ) return HX;
        if(second == HY) return H;
    }

    if(first == Z){
    
        if(second == X) return Y;
        if(second == Y) return X;
        if(second == Z) return I;
        if(second == H) return HX;
        if(second == HX) return H;
        if(second == HZ) return HY;
        if(second == HY) return HZ;
    
    }
    if(first == H){
        if(second == X) return HX;
        if(second == Y) return HY;
        if(second == Z) return HZ;
        if(second == H) return I;
        if(second == HX) return X;
        if(second == HZ) return Z;
        if(second == HY) return Y;
    }

    if(first == HX){
        if(second == X) return H;
        if(second == Y) return HZ;
        if(second == Z) return HY;
        if(second == H) return Z;
        if(second == HX) return Y;
        if(second == HZ) return I;
        if(second == HY) return X;
    }

    if(first == HY){
    
        if(second == X) return HZ;
        if(second == Y) return H;
        if(second == Z) return HX;
        if(second == H) return Y;
        if(second == HX) return Z;
        if(second == HZ) return X;
        if(second == HY) return I;
    }

    if(first == HZ){
    
        if(second == X) return HY;
        if(second == Y) return HX;
        if(second == Z) return H;
        if(second == H) return X;
        if(second == HX) return I;
        if(second == HZ) return Y;
        if(second == HY) return Z;
    }
}


