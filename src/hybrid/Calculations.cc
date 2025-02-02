#include <iostream>
#include <iomanip>  // this include file enables the use of I/O manipulators for controlling the format of the output stream
#include <fstream>
#include <complex>
#include <cstdarg>  // this include file enables working with variable length parameter lists

using namespace std;

#include "diracMatrices.h"
#include "Matrix.h"
#include "Calculations.h"
#include "Constants.h"
#include "Operators.h"
#include "ABCDE.h"

#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */




double Do_Fivefold_Calc( Reaction_parameters *Reac, complex<double> Lepton_T[][4], Hadron_prime * Had)
{

// // // // NOTATION // // // // 
// process = 1   CC interaction
//   HELICITY == -1 --> W^+ induced 1-pion production 
//  nucleon = 1 --> proton initial state
//    decay = 1 or 2 --> p + pi^+
//  nucleon = 2 --> neutron initial state
//    decay = 1 --> p + pi^0
//    decay = 2 --> n + pi^+
//   
//   HELICITY == 1 --> W^- induced 1-pion production 
//  nucleon = 2 --> neutron initial state
//    decay = 1 or 2 --> n + pi^-
//  nucleon = 1 --> proton initial state
//    decay = 1 --> n + pi^0
//    decay = 2 --> p + pi^+
// //    
// process = 0   EM interaction --> photon induced 1-pion production 
//  nucleon = 1 --> proton initial state 
//    decay = 1 --> p + pi^0
//    decay = 2 --> n + pi^+
//  nucleon = 2 --> neutron initial state
//    decay = 1 --> n + pi^0
//    decay = 2 --> p + pi^-
// // // // // // // // // // //


  /******************************************
   *       Construct the lepton tensor      *
   ******************************************/
  // In this version the Lepton tensor is already computed and stored in the struct 'Lep'

  /******************************************
   *       Construct the Operators          *
   ******************************************/

  Matrix Op[4];
  GetOperator(Had, Reac, Op); // The type of operator depends on the flag iModel which is set in constants.h!


  
  /***************************
   * BUILD THE HADRON TENSOR *
   ***************************/

  // Hadron tensor is computed by using the Dirac identity for free on-shell spinors

  // Building the matrices (kslash + M)*Gamma
  
  complex<double> block1[4][4], block2[4][4];
  
  double Mfac;
  for (int i = 0 ; i < 4 ; i++)
  {
    for (int j = 0 ; j < 4 ; j++)
    {
      Mfac = real(MN*Id.M[i][j]);
      block1[i][j] = (1./(4.*MN))*(Had->kiSlash.M[i][j] + Mfac)*Gamma[0].M[j][j];
      block2[i][j] = (1./(2.*MN))*(Had->kNSlash.M[i][j] + Mfac)*Gamma[0].M[i][i];
    }
  }
  // Block 1 and 2 contain the factors 1/2 and 1/4 from summation over helicities already, they are not added later.

  // Conjugating the Operator, and then computing (ki_slash + M)*Op^+*(kN_slash +M)
  Matrix Conj_Op[4];
  for (int i = 0 ; i < 4 ; i++)
  {
    Conj_Op[i] = ConjMatrix(Op[i]);
  }

  Matrix Sandwich_Op[4];
  for (int i = 0 ; i < 4 ; i++)
  {
    for (int m = 0 ; m < 4 ; m++)
    {
      for (int j = 0 ; j < 4 ; j++)
      {
        Sandwich_Op[i].M[j][m] = 0.0; 
        for (int k = 0 ; k < 4 ; k++)
        {
          for (int l = 0 ; l < 4 ; l++)
          {
            Sandwich_Op[i].M[j][m] += block1[j][l]*Conj_Op[i].M[l][k]*block2[k][m];
          }
        }
        
      }
    }
  }

  complex<double> Hadron[4][4] = {}; // Hadron tensor is Trace(Sandwich[i]*Op[j])

  for(int i=0; i<4; i++)
  {
    for(int j=0; j<4; j++)
    {
      Trace_product(Sandwich_Op[i], Op[j], Hadron[i][j]); // The trace op Sandwich_Op*Op is stored in Hadron[i][j]
    }
  }

  // Calculate the factors ABCDE, without any prefactors except for the additional normalization in Lepton tensor and 1/(8MN^2) from hadron
  ABCDE(Hadron, Lepton_T, Reac->ABCDE, Reac->R_factors);

  return 0;
}


