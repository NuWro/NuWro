#ifndef DELTA_H
#define DELTA_H
#include "Matrix.h"
#include <complex>

// // Oset and Salcedo medium modification of the delta width
void OSMM( int medmod, double s, double r, complex<double> &f_OSMM);

void Delta_ff(double s, double u, double &Deltaff);


void Gamma_WNDelta(int cross, int process,  double Qsq, double s, double Q[], double kResonance[], double ki[], Matrix WNDelta[][4] );

/***********************************************************************
                         THE DELTA PROPAGATOR
***********************************************************************/
void S_Delta( int medmod, int cross, double W2, double kResonance[], Matrix kRSlash, Matrix Delta[][4] );

void Gamma_DeltaNpi( double s, double fDeltaNpi, double kResonance[], double kpi[], Matrix DeltaNpi[] );
void DP_current( int medmod, int Pascalutsa, int process, int nucleon, int decay, int Helicity, int cross, double Qsq, double s, double u, double Q[], double kResonance[], double ki[], double kpi[], Matrix kRSlash, Matrix Op_delta[] );

void DP_current_pre( int medmod, int Pascalutsa, int process, int nucleon, int decay, int Helicity, int cross, double Qsq, double s, double u, double Q[], double kResonance[], double ki[], double kpi[], Matrix kRSlash,
 const Matrix WNDelta[][4],
 const Matrix PropDelta[][4],
Matrix Op_delta[] );

void CDP_current( int medmod, int Pascalutsa, int process, int nucleon, int decay, int Helicity, int cross, double Qsq, double s, double u, double minusQ[], double kResonance[], double kN[], double kpi[], Matrix kRSlash, Matrix Op_delta_cross[] );

#endif
