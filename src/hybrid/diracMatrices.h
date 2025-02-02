#ifndef DIRAC_MATRICES
#define DIRAC_MATRICES

#include "Matrix.h"
/** 
 * this is really not very clean, introducing global variables and stuff
 * but raul wants it that way.
 */
extern const Matrix Id;
extern const Matrix Gamma[4];
extern const Matrix Gamma5;
extern const Matrix Gamma_mu5[4];
extern const Matrix Gamma_munu[4][4];
extern const Matrix mGamma_munu[4][4];




#endif // DIRAC_MATRICES
