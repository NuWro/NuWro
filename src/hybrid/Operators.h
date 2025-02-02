#ifndef OPERATORS_H
#define OPERATORS_H
#include "Constants.h"
#include "Matrix.h"

void GetOperator(Hadron_prime *kin, Reaction_parameters *par, Matrix Operator[4]);

void Get_Res_ChPT(Hadron_prime *kin, Reaction_parameters *params, Matrix ChPT[4], Matrix Resonances[4]);

void Get_ReChi(Hadron_prime *kin, Reaction_parameters *params, Matrix ReChi[4]);


#endif
