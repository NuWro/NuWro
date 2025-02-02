#ifndef CONSTANTS_H
#define CONSTANTS_H


/*
The "Constants.h" header file, when included, provides some useful constant objects that will be used throughout the code.
*/

#include <complex>
#include "Matrix.h"
using namespace std;

// The complex number unit
const complex<double> I = complex<double> (0, 1);

// Some masses
const double Mp = 938.27203;  // proton mass (MeV)
const double Mn = 939.56536;  // neutron mass (MeV)
const double MN = 938.918695;  // average mass (MeV)
const double MN2 = 881568.315821;//MN^2

const double Mpi_chrgd = 139.57018;  // charged-pion mass (MeV)
const double Mpi_ntrl = 134.9766;  // neutral pion mass (MeV)
const double Mpi = 138.0389867; // average (2*Mpi_chrgd+Mpi_ntrl)/3
const double Mpi2 = 19054.761849; //Mpi^2

const double Pi = 3.141592654;  // \pi
const double hbc = 197.3270;  // hbar*c (MeV*fm)


// Electroweak constants
const double M_W = 80385.;  // W-boson mass (MeV)
const double M_Z = 91187.6; // Z-boson mass (MeV)
const double Cabibbo = 0.974;  // cosine Cabibbo angle
const double G_Fermi = 1.16637e-11;  // Fermi constant (MeV^{-2})
const double Alpha = 0.007297353;  // fine structure constant
const double Sin2W = 0.23122; 
const double QWeak = 0.07512; //weak charge of the proton = 1-4*Sin2W

// Delta-related constants
const double Wdth = 120.;  	// Delta-particle width (MeV)
const double MDelta = 1232.;  	// Delta-particle mass (MeV)
const double fDelta = 2.18; 	// 

// D13-related constants 
const double WdthD13 = 115.; 
const double MD13 = 1515.; 
const double fD13 = 1.62; // br = 0.60

// // S11-related constants 
const double WdthS11 = 150.; 
const double MS11 = 1535.; 
const double fS11 = 0.16; // br = 0.45

// // P11-related constants 
const double WdthP11 = 350.;
const double MP11 = 1430.; 
const double fP11 = 0.489; // br = 0.60 

// // Non-resonant constants 
const double gA = 1.26;
const double fpi = 93.;

///////////////////////////////////////////////////////////
// //    PARAMETERS RELATED WITH THE OPERATOR         // //
///////////////////////////////////////////////////////////
//Add as an input later!
//const int iModel = 0; // iModel=0 -> LEM without form factors for the resonances, iModel=1 -> LEM with resonance form factors removed; iModel=2 (regge) removed

//Parameters used for the operator which are not kinematic;
//const int Pascalutsa = 0;// 1 to use Pascalutsa Delta decay vertex, '0' for the normal one (Delta pole)
//int medmod = 0; //Medmod is in reaction parameters
// // // // // // // // // // // // // // // // // // // // // // // // // 

// for ease of writing
const int down = 0; //spin -1/2
const int up = 1;   //spin 1/2

const int proton = 1;
const int neutron = 2;

// Structs to pass information around

struct Reaction_parameters
{
  // parameters for type of reaction/Operator
  int process, decay, Helicity, medmod, nucleon; 
  int Pauli_blocking;
  int irun;
  double Minf_A, lambda_rho;

  // Output
  double *ABCDE, *R_factors;
};


struct Hadron_prime
{
  double *Q, *minusQ;
  double *kpi, *kN, *ki;         // 4vectors in the kN//z system

  Matrix kpiSlash, QSlash, kiSlash, kNSlash; 

  double *sMan, *tMan, *uMan;
  Matrix sSlash, tSlash, uSlash;

  double W, s, t, u, Qsq;                     //invariants

  double k_pi;  // Scalars
};


#endif
