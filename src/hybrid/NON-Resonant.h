#ifndef NONRESONANT_H
#define NONRESONANT_H
#include "Constants.h"
#include "Matrix.h"
#include "diracMatrices.h"
// /***************************************************************************
//                           THE -W N N- VERTEX
// ***************************************************************************/

void Gamma_WNN( int nucleon, int process, int decay, int cross, double Qsq, double Q[], const Matrix& Qslash, Matrix WNN[] );


/***********************************************************************
                         THE Nucleon PROPAGATOR
***********************************************************************/


void S_Nprop( double W2, const Matrix& kSlash, Matrix &Nprop );


void NP_current( int process, int nucleon, int decay, int Helicity, int cross, double Qsq, double W2, double Q[], const Matrix& Qslash, const Matrix& sSlash, const Matrix& kpiSlash, Matrix Op_NP[] );


void CNP_current( int process, int nucleon, int decay, int Helicity, int cross, double Qsq, double W_cross2, double Q[], const Matrix& Qslash, const Matrix& uSlash, const Matrix& kpiSlash, Matrix Op_NP_cross[] );


#endif
