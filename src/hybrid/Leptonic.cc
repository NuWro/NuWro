/*****************************************************************************

*                            THE LEPTON TENSOR                               *

*****************************************************************************/

#include <stdlib.h>
#include <math.h>

void Leptonic_Tensor_fourvector(const int &process, double kl[4], double kl_inc[4], 
		      double Lepton_S[][4], double Lepton_A[][4] )
{
  //Identical to the original function, but takes the four vectors as arguments and does not include the energy normalization that cancels with prefactors anyway, -> Transforms correct under LT

  double kl_inc_dot_kl = kl_inc[0]*kl[0] - kl_inc[1]*kl[1] - kl_inc[2]*kl[2] - kl_inc[3]*kl[3];

  double Norm_factor;
  
  if( process == 0 ){
    Norm_factor = 1./4;  
  }else{
    Norm_factor = 1.;  
  }

  //NOTE: No energy dependence in the norm factor, to restore lorentz invariance of matrix element defined as Leptonic*Hadron contraction in Calculations.cpp, 

  
  Lepton_S[0][0] = Norm_factor*(2.*kl_inc[0]*kl[0] - kl_inc_dot_kl);
  Lepton_S[0][1] = Norm_factor*(-kl_inc[0]*kl[1] + -kl[0]*kl_inc[1]);  // mind the sign difference for the spatial compts!
  Lepton_S[0][2] = Norm_factor*(-kl_inc[0]*kl[2] + -kl[0]*kl_inc[2]);
  Lepton_S[0][3] = Norm_factor*(-kl_inc[0]*kl[3] + -kl[0]*kl_inc[3]);

  Lepton_S[1][0] = Lepton_S[0][1];
  Lepton_S[1][1] = Norm_factor*(2*kl_inc[1]*kl[1] + kl_inc_dot_kl);
  Lepton_S[1][2] = Norm_factor*(kl_inc[1]*kl[2] + kl[1]*kl_inc[2]);
  Lepton_S[1][3] = Norm_factor*(kl_inc[1]*kl[3] + kl[1]*kl_inc[3]);

  Lepton_S[2][0] = Lepton_S[0][2];
  Lepton_S[2][1] = Lepton_S[1][2];
  Lepton_S[2][2] = Norm_factor*(2*kl_inc[2]*kl[2] + kl_inc_dot_kl);
  Lepton_S[2][3] = Norm_factor*(kl_inc[2]*kl[3] + kl[2]*kl_inc[3]);

  Lepton_S[3][0] = Lepton_S[0][3];
  Lepton_S[3][1] = Lepton_S[1][3];
  Lepton_S[3][2] = Lepton_S[2][3];
  Lepton_S[3][3] = Norm_factor*(2*kl_inc[3]*kl[3] + kl_inc_dot_kl);

  
    Lepton_A[0][0] = 0.;
    Lepton_A[0][1] = Norm_factor*(kl_inc[2]*kl[3] - kl[2]*kl_inc[3]);
    Lepton_A[0][2] = Norm_factor*(-kl_inc[1]*kl[3] + kl[1]*kl_inc[3]);
    Lepton_A[0][3] = Norm_factor*(kl_inc[1]*kl[2] - kl[1]*kl_inc[2]);
    
    Lepton_A[1][0] = -Lepton_A[0][1];
    Lepton_A[1][1] = 0.;
    Lepton_A[1][2] = Norm_factor*(kl_inc[0]*kl[3] - kl[0]*kl_inc[3]);
    Lepton_A[1][3] = Norm_factor*(-kl_inc[0]*kl[2] + kl[0]*kl_inc[2]);
    
    Lepton_A[2][0] = -Lepton_A[0][2];
    Lepton_A[2][1] = -Lepton_A[1][2];
    Lepton_A[2][2] = 0.;
    Lepton_A[2][3] = Norm_factor*(kl_inc[0]*kl[1] - kl[0]*kl_inc[1]);
    
    Lepton_A[3][0] = -Lepton_A[0][3];
    Lepton_A[3][1] = -Lepton_A[1][3];
    Lepton_A[3][2] = -Lepton_A[2][3];
    Lepton_A[3][3] = 0.;


  
}

