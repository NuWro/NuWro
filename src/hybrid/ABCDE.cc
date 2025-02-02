#include <complex>
using namespace std;
#include "Constants.h"

void ABCDE(const std::complex<double> Hadron[4][4], const std::complex<double> Lepton[4][4], double *ABC_factors, double *R_factors)
{
  // We calculate the factors A B C D and E from lepton and hadron tensor as described in eq. 31 from Sobscyk paper.
  // Note a different normalization of lepton and hadron tensors in our convention
  // We store them in the array ABC factors
  double Hadron_s[4][4];
  double Hadron_a[4][4]; // Symmetric and antisymmetric hadron tensor elements. 
  for (int mu = 0 ; mu < 4 ; mu++)
  {
    for (int nu =0 ; nu < 4 ; nu++)
    {
      Hadron_s[mu][nu] = real(0.5*(Hadron[mu][nu]+Hadron[nu][mu]));
      Hadron_a[mu][nu] = real(-I*0.5*(Hadron[mu][nu]-Hadron[nu][mu]));
    }
  }


  // A:
  ABC_factors[0] = real(Lepton[0][0]*Hadron_s[0][0] + 2.*Lepton[3][0]*Hadron_s[3][0] + Lepton[3][3]*Hadron_s[3][3]
                 + 0.5*(Lepton[1][1]+Lepton[2][2])*(Hadron_s[1][1]+Hadron_s[2][2]) + 2.*I*Lepton[1][2]*Hadron_a[1][2]);

  // B:
  ABC_factors[1] = real( 2.*(Lepton[0][1]*Hadron_s[0][1]+Lepton[1][3]*Hadron_s[1][3]  + I*Lepton[0][2]*Hadron_a[0][2] + I*Lepton[2][3]*Hadron_a[2][3]));

  // C:
  ABC_factors[2] = real(0.5*(Lepton[1][1]-Lepton[2][2])*(Hadron_s[1][1]-Hadron_s[2][2]));

  // D:
  ABC_factors[3] = real(2.*(-Lepton[0][1]*Hadron_s[0][2] -Lepton[1][3]*Hadron_s[2][3] + I*Lepton[0][2]*Hadron_a[0][1]+ I*Lepton[2][3]*Hadron_a[1][3]));

  // E:
  ABC_factors[4] = real((Lepton[2][2] - Lepton[1][1])*Hadron_s[1][2]);

  // Terms for calculation of the inclusive cross section (A structure function)
  R_factors[0] = Hadron_s[0][0];
  R_factors[1] = Hadron_s[3][0];
  R_factors[2] = Hadron_s[3][3];
  R_factors[3] = Hadron_s[1][1]+Hadron_s[2][2];
  R_factors[4] = Hadron_a[1][2];

}


