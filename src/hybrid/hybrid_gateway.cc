#include <array>
#include <complex>
#include <cstdarg>

using namespace std;

#include "Matrix.h"
#include "Calculations.h"
#include "Constants.h"
#include "diracMatrices.h"
#include "Leptonic.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
// // // // // // // // // // // // // // 

#include <cstdlib>
#include <cassert>
#include <cstdio>
#include <fstream>

// int invcubic_real(double a, double b, double c, double d, double& sol) 
// {
//     //Solving a cubic equation using the discriminant approach
//     //We only need real solutions, no complex parts are calculated.
//     //Moreover, the real solution should be between -1 and 1 (this is not garanteed!)
//     b /= a;
//     c /= a;
//     d /= a;
    
//     double disc, q, r, dum1, s, t, term1, r13;
//     q = (3.0*c - (b*b))/9.0;
//     r = -(27.0*d) + b*(9.0*c - 2.0*(b*b));
//     r /= 54.0;
//     disc = q*q*q + r*r;
//     term1 = (b/3.0);
    
//     double x1_real, x2_real, x3_real;
//     if (disc > 0)   // One root real, two are complex
//     {
//         s = r + sqrt(disc);
//         s = s<0 ? -cbrt(-s) : cbrt(s);
//         t = r - sqrt(disc);
//         t = t<0 ? -cbrt(-t) : cbrt(t);
//         sol = -term1 + s + t;
//     } 
//     // The remaining options are all real
//     else if (disc == 0)  // All roots real, at least two are equal.
//     { 
//         r13 = r<0 ? -cbrt(-r) : cbrt(r);
//         x1_real = -term1 + 2.0*r13;
//         if (abs(x1_real) <= 1)
//             {
//                 sol = x1_real;

//             }else{
//                 sol = -(r13 + term1);
//             }
//     }
//     // Only option left is that all roots are real and unequal (to get here, q < 0)
//     else
//     {
//         q = -q;
//         dum1 = q*q*q;
//         dum1 = acos(r/sqrt(dum1));
//         r13 = 2.0*sqrt(q);
//         x1_real = -term1 + r13*cos(dum1/3.0);
//         x2_real = -term1 + r13*cos((dum1 + 2.0*Pi)/3.0);
//         if (abs(x1_real) <= 1.){sol = x1_real;}
//         else if( abs(x2_real) <= 1){sol = x2_real;}
//         else{
//            sol = -term1 + r13*cos((dum1 + 4.0*Pi)/3.0);
//         }

//     }
//     return 0;  
// }

// void Calc_parabola(double* x, double* y, double &A, double &B, double &C)
// {
//     //Calculate the constants in A*x^2 + B*x + C = 0 going through 3 points whose x,y values are in the input arrays

//     double denom = (x[0] - x[1]) * (x[0] - x[2]) * (x[1] - x[2]);
//     A     = (x[2] * (y[1] - y[0]) + x[1] * (y[0] - y[2]) + x[0] * (y[2] - y[1])) / denom;
//     B     = (x[2]*x[2] * (y[0] - y[1]) + x[1]*x[1] * (y[2] - y[0]) + x[0]*x[0] * (y[1] - y[2])) / denom;
//     C     = (x[1] * x[2] * (x[1] - x[2]) * y[0] + x[2] * x[0] * (x[2] - x[0]) * y[1] + x[0] * x[1] * (x[0] - x[1]) * y[2]) / denom;

// }

std::array<double, 5> get_lepton_vec(double El_inc, double Q2, double W,
                                     double leptonmass, double nucleonmass,
                                     int *params) {
  std::array<double, 5> lepton_vec;
  // Effective nucleon mass
  double nucleonmass2 = nucleonmass * nucleonmass;

  // Making the struct, could be removed at some point, because it may not be
  // efficient
  Hadron_prime Hadronbag{};

  // Set arrays for lepton kinematics
  double kl[4]{}, kl_inc[4]{}, Q[4]{}, Pa[4]{};
  complex<double> Lepton_T[4][4]{}; // lepton tensor

  Reaction_parameters Reacbag{};
  double ABCDE[5]{}, R_factors[5]{};
  // Fill Reac
  Reacbag.process = params[0], Reacbag.decay = params[1],
  Reacbag.Helicity = params[2], Reacbag.nucleon = params[3];
  Reacbag.ABCDE = ABCDE, Reacbag.R_factors = R_factors;

  // Getting the pointers:
  //  Hadron_prime* Had = &Hadronbag;
  Reaction_parameters *Reac = &Reacbag;

  //////////////////////////////////////////////
  ///////Lepton kinematics:
  //////////////////////////////////////////////////

  double W2 = W * W;
  double QQ = Q2;

  double Enu = El_inc;
  double w =
      (W * W + QQ - nucleonmass2) / (2 * nucleonmass); // This omega is  in LAB!
  // if (w <= 0 || w > Enu){return -1;}
  // if (pow((Enu - w),2) - leptonmass*leptonmass <= 0){return -1;}
  double k_l = sqrt(pow((Enu - w), 2) - leptonmass * leptonmass);
  double E_l = Enu - w;
  // if (fabs((QQ+w*w - pow(Enu,2) - pow(k_l,2))/(2.*Enu*k_l)) > 1){return -1;}
  double radThetal =
      acos((QQ + w * w - pow(Enu, 2) - pow(k_l, 2)) / (-2. * Enu * k_l));
  double q = sqrt(QQ + w * w);

  // define the boost parameters
  double v = q / (w + nucleonmass);
  double gamma = (w + nucleonmass) / W;

  // define CMS angles for the leptons
  double costheta = (Enu - k_l * cos(radThetal)) / q;
  double sintheta = k_l * sin(radThetal) / q;

  // The leptonic side in CMS, q~z and in x-z plane is now given as:

  kl[0] = gamma * (E_l - v * (Enu * costheta - q));
  kl[1] = sintheta * Enu;
  kl[2] = 0.;
  kl[3] = gamma * (-1 * v * E_l + (Enu * costheta - q));

  kl_inc[0] = gamma * (El_inc - v * Enu * costheta);
  kl_inc[1] = kl[1];
  kl_inc[2] = 0;
  kl_inc[3] = gamma * (-1 * v * Enu + Enu * costheta);

  Q[0] = kl_inc[0] - kl[0];
  Q[1] = 0;
  Q[2] = 0;
  Q[3] = kl_inc[3] - kl[3];

  double Lepton_S[4][4]{}, Lepton_A[4][4]{};
  Leptonic_Tensor_fourvector(Reac->process, kl, kl_inc, Lepton_S, Lepton_A);
  double h = Reac->Helicity;
  double CC = Reac->process;
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      Lepton_T[i][j] = Lepton_S[i][j] - I * Lepton_A[i][j] * h * CC;
    }
  }
  lepton_vec[0] = real(Lepton_T[0][0]);
  lepton_vec[1] = real(Lepton_T[3][0]);
  lepton_vec[2] = real(Lepton_T[3][3]);
  lepton_vec[3] = real(Lepton_T[1][1] + Lepton_T[2][2]);
  lepton_vec[4] = -real(I * Lepton_T[1][2]);
  return lepton_vec;
}

int hybrid_ABCDE(double El_inc, double Q2, double W, double leptonmass, double nucleonmass, double *costheta_pi_in, int N, int *params, double (*strucfuncs)[5], double (*Inclusive)[5])
{
  // Effective nucleon mass
  double nucleonmass2 = nucleonmass*nucleonmass;

  //Making the struct, could be removed at some point, because it may not be efficient
  Hadron_prime Hadronbag{};

  //Set arrays for lepton kinematics
  double kl[4]{},kl_inc[4]{}, Q[4]{}, Pa[4]{};
  complex<double> Lepton_T[4][4]{}; // lepton tensor

  //Same for Hadron Kin
  double xQ[4]{}, xminusQ[4]{};
  double kpi[4]{}, kN[4]{}, ki[4]{};         //4vectors in the kN//z system
  Hadronbag.Q= xQ, Hadronbag.minusQ=xminusQ, Hadronbag.kpi = kpi, Hadronbag.kN=kN, Hadronbag.ki=ki;

  //Mandelstam variables
  double sMan[4]{}, tMan[4]{}, uMan[4]{};
  Hadronbag.sMan=sMan;
  Hadronbag.tMan=tMan;
  Hadronbag.uMan=uMan;

  Reaction_parameters Reacbag{};
  double ABCDE[5]{}, R_factors[5]{};
  //Fill Reac
  Reacbag.process = params[0], Reacbag.decay= params[1], Reacbag.Helicity = params[2] , Reacbag.nucleon = params[3];
  Reacbag.ABCDE = ABCDE, Reacbag.R_factors = R_factors;

  //Getting the pointers:
  Hadron_prime* Had = &Hadronbag;
  Reaction_parameters* Reac = &Reacbag;

  //////////////////////////////////////////////
  ///////Lepton kinematics:
  //////////////////////////////////////////////////

  double W2 = W*W;
  double QQ = Q2;

  double Enu = El_inc;
  double w = (W*W + QQ - nucleonmass2)/(2*nucleonmass); //This omega is  in LAB!
  if (w <= 0 || w > Enu){return -1;}
  if (pow((Enu - w),2) - leptonmass*leptonmass <= 0){return -1;}
  double k_l = sqrt(pow((Enu - w),2) - leptonmass*leptonmass);
  double E_l = Enu - w;
  if (fabs((QQ+w*w - pow(Enu,2) - pow(k_l,2))/(2.*Enu*k_l)) > 1){return -1;}
  double radThetal = acos((QQ+w*w - pow(Enu,2) - pow(k_l,2))/(-2.*Enu*k_l) );
  double q = sqrt(QQ + w*w);

  //define the boost parameters
  double v = q/(w+nucleonmass);
  double gamma = (w + nucleonmass)/W;

  //define CMS angles for the leptons
  double costheta = (Enu - k_l*cos(radThetal))/q;
  double sintheta = k_l*sin(radThetal)/q;

  //The leptonic side in CMS, q~z and in x-z plane is now given as:

  kl[0] = gamma*(E_l - v*(Enu*costheta - q));
  kl[1] = sintheta*Enu;
  kl[2] = 0.;
  kl[3] = gamma*(-1*v*E_l + (Enu*costheta - q));

  kl_inc[0] = gamma*(El_inc - v*Enu*costheta); 
  kl_inc[1] = kl[1];
  kl_inc[2] = 0;
  kl_inc[3] = gamma*(-1*v*Enu+ Enu*costheta);

  Q[0] = kl_inc[0] - kl[0];
  Q[1] = 0;
  Q[2] = 0;
  Q[3] = kl_inc[3] - kl[3];

  double Lepton_S[4][4]{}, Lepton_A[4][4]{};
  Leptonic_Tensor_fourvector(Reac->process,kl, kl_inc, Lepton_S, Lepton_A);
  double h = Reac->Helicity;
  double CC = Reac->process;
  for (int i = 0 ; i < 4 ; i++){
    for (int j = 0; j < 4 ; j++)
        {
      Lepton_T[i][j] = Lepton_S[i][j] - I*Lepton_A[i][j]*h*CC;
        }
  }
  /////////////////////////////////////////////////////
  ///////HADRON kinematics (the part that does not depend on cosine of pion)
  ////////////////////////////////////////////////////////////

  double E_pi = (W*W + pow(Mpi,2) - nucleonmass2)/(2*W);
  double E_N = W - E_pi;
  double p_N = nucleonmass*q/W;

  Had->Q[0] = Q[0];
  Had->Q[1] = Q[1];
  Had->Q[2] = Q[2];
  Had->Q[3] = Q[3];


  Had->minusQ[0] = -Had->Q[0];
  Had->minusQ[1] = -Had->Q[1];
  Had->minusQ[2] = -Had->Q[2];
  Had->minusQ[3] = -Had->Q[3];

  //the pion and final nucleon 
  double Tpi = E_pi - Mpi;
  if (Tpi < 0){return -1;}

  Had->ki[0] = W - Had->Q[0];
  Had->ki[1] = 0;
  Had->ki[2] = 0;
  Had->ki[3] = -p_N;

  Had->kiSlash =  Had->ki[0]*Gamma[0]-Had->ki[3]*Gamma[3]; //Usin explicitly that x-y components of ki are zero
  Had->QSlash = Had->Q[0]*Gamma[0]-Had->Q[3]*Gamma[3];

  //Mandelstam s is trivial in CMS:
  Had->sMan[0] = W;
  Had->sMan[1] = 0;
  Had->sMan[2] = 0;
  Had->sMan[3] = 0;
  Had->sSlash = Had->sMan[0]*Gamma[0]; 

  Had->Qsq = QQ;

  Had->k_pi = sqrt(Tpi*(Tpi+2*Mpi));
  Had->kpi[0] = E_pi;
  Had->kpi[2] = 0; //pion in x-z frame

  Had->kN[0] = W - E_pi;
  if (Had->kN[0] < nucleonmass){return -1;}
  Had->kN[2] = 0;

  Had->W = W;
  Had->s = W*W;

  for (int i = 0 ; i < N ; i++)
  {
    double CT = costheta_pi_in[i];
    double radThetapi = acos(CT);

    Had->kpi[1] = sin(radThetapi)*Had->k_pi;
    Had->kpi[3] = cos(radThetapi)*Had->k_pi;

    Had->kN[1] = -Had->kpi[1];
    Had->kN[3] = -Had->kpi[3];

    double prefactor = 0.5;

    //Mandelstam
    for(int i=0; i<4; i++)
    {
      Had->uMan[i] = Had->kN[i]  - Had->Q[i];
      Had->tMan[i] = Had->Q[i]   - Had->kpi[i];
    }
    Had->u = pow(Had->uMan[0],2) - pow(Had->uMan[1],2) - pow(Had->uMan[2],2) - pow(Had->uMan[3],2);
    Had->t = pow(Had->tMan[0],2) - pow(Had->tMan[1],2) - pow(Had->tMan[2],2) - pow(Had->tMan[3],2);

    //Slashes (to be optimized!)
    Had->kNSlash = Had->kN[0]*Gamma[0]-Had->kN[1]*Gamma[1]-Had->kN[2]*Gamma[2]-Had->kN[3]*Gamma[3];
    Had->kpiSlash = Had->kpi[0]*Gamma[0]-Had->kpi[1]*Gamma[1]-Had->kpi[2]*Gamma[2]-Had->kpi[3]*Gamma[3];
    Had->tSlash = Had->tMan[0]*Gamma[0] - Had->tMan[1]*Gamma[1]- Had->tMan[2]*Gamma[2] - Had->tMan[3]*Gamma[3];
    Had->uSlash = Had->uMan[0]*Gamma[0] - Had->uMan[1]*Gamma[1]- Had->uMan[2]*Gamma[2] - Had->uMan[3]*Gamma[3];

    Do_Fivefold_Calc(Reac, Lepton_T, Had);
    if(Inclusive)
    {
      Inclusive[i][0] = Reac->R_factors[0];
      Inclusive[i][1] = Reac->R_factors[1];
      Inclusive[i][2] = Reac->R_factors[2];
      Inclusive[i][3] = Reac->R_factors[3];
      Inclusive[i][4] = Reac->R_factors[4];
    }
    if(strucfuncs)
    {
      strucfuncs[i][0] = Reac->ABCDE[0] * prefactor;
      strucfuncs[i][1] = Reac->ABCDE[1] * prefactor;
      strucfuncs[i][2] = Reac->ABCDE[2] * prefactor;
      strucfuncs[i][3] = Reac->ABCDE[3] * prefactor;
      strucfuncs[i][4] = Reac->ABCDE[4] * prefactor;
    }



  }

  return 0;
}

// int genevents(double *x, double *y, int Nevents, double *events, double &A, double &B, double &C)
// {
// //Calculate N values of cos\\theta according to the distributions given by a parabola going through 3 points determined by x and y values in the arrays x,y by inversion of CDF
//     //(output) A, B, C coefficients for the PDF: A*x**2 + B*x + C
//     double a,b,c,d; //coefficients for the CDF : a*x**3 + b*x**2 + c*x + d, normalized to one in the interval [-1:1]
//     double int_pdf; //value of the integral of PDF from [-1:1]

//     //Calculate the constants in A*x^2 + B*x + C = 0 going through 3 points whose x,y values are in the input arrays
//     double denom = (x[0] - x[1]) * (x[0] - x[2]) * (x[1] - x[2]);
//     A     = (x[2] * (y[1] - y[0]) + x[1] * (y[0] - y[2]) + x[0] * (y[2] - y[1])) / denom;
//     B     = (x[2]*x[2] * (y[0] - y[1]) + x[1]*x[1] * (y[2] - y[0]) + x[0]*x[0] * (y[1] - y[2])) / denom;
//     C     = (x[1] * x[2] * (x[1] - x[2]) * y[0] + x[2] * x[0] * (x[2] - x[0]) * y[1] + x[0] * x[1] * (x[0] - x[1]) * y[2]) / denom;

//     //Integral of the parabola over [-1:1] is:
//     int_pdf = 2.*A/3. + 2.*C;

//     // The a,b,c,d for the CDF : cubic(x) determined by the integral over the parabola from [-1:x]
//     a = A/3.;
//     b = B/2.;
//     c = C;
//     d = -1.*(-1.*A/3. + B/2. -C);

//     //normalize these such that the CDF is one at x=1
//     a/=int_pdf;
//     b/=int_pdf;
//     c/=int_pdf;
//     d/=int_pdf;

//     //We invert the CDF for Nevents random numbers, rand, between 0 and 1
//     for (int i=0; i < Nevents; i++)
//     {
//         double sol; //stores the solution
//         double r = ( (double) rand()/ (RAND_MAX)); //random number between 0 and 1

//         //Compute the real solution of
//         //a*x**3 + b*x**2 + c*x + (d-r) = 0
//         //The solution that lies between -1 and 1 is stored in sol
//         invcubic_real(a, b, c, d-r,sol); 

//         if (abs(sol) <= 1)
//         {
//             events[i]=sol;
//         }else{
//             //Sol is not between 0 and 1!, this should never happen, It can only happen due to possibly very small numerical inacurracy in the cubic roots.
//             //If it happens, I just ignore this value and make an extra event
//             cout << "Invalid value in inversion! = " << sol << endl;
//             i--;
//         }

//     }
//     return 0;
// }

// int main()
// {
// double leptonmass = 0.0;
// ofstream histogram; // an outputfile for the raw events of costhete
// ofstream distribution; // an outputfile for the theoretical distributions etc

// int params[4];
// params[0] = 1; //process
// params[1]=2; //decay
// params[2]=-1; //Helicity
// params[3] = 2; //nucleon
// double El_inc = 5000.0;
// double Q2 = 0.2*1e6;
// double W = 1230;
// double costheta_pi_in[2000];
// int cnt = 0;
// for (double i = -1. ; i <= 1. ; i+=0.01)
// {
// 	costheta_pi_in[cnt] = i;	
// 	cnt++;
// }
// int N = cnt;
// double strucfuncs[cnt][5];
// double Inclusive[cnt][5];

// double cospi_sample[3];
// double strucfuncs_sample[3][5];
// cospi_sample[0]=-0.75;
// cospi_sample[1]=0.;
// cospi_sample[2]=0.75;

// double E_pi = (W*W + pow(Mpi,2) - pow(MN,2))/(2*W);
// double kpi = sqrt(E_pi*E_pi - Mpi*Mpi);

// double F = pow(G_Fermi*Cabibbo,2)/2.;
// double A,B,C;
// double y[3];

// for (double Q2= 0.1*1e6 ; Q2 < 3.0*1e6 ; Q2 += 0.5*1e6)
// {
//     for (double W = 1100 ; W < 1500 ; W +=66)
//     {
//         //output of the events in file:
// //        std::string histfnm = "hist_"+to_string(params[3]) + to_string(params[1]) + "_" + to_string((int)W) + "_" + to_string((int)(Q2*1e-3)) +".dat"; //nasty cpp

//         //output of the functions in file:
//         std::string funcfnm = "cos_"+to_string(params[3]) + to_string(params[1]) + "_" + to_string((int)W) + "_" + to_string((int)(Q2*1e-3)) +".dat"; //nasty cpp
        
//         double integrated = 0;
//         double chi2 = 0;
//         int ERR =  hybrid_ABCDE(El_inc, Q2, W, leptonmass ,costheta_pi_in, N, params, strucfuncs,Inclusive);
//         if (ERR == -1){continue;}
//         ERR =  hybrid_ABCDE(El_inc, Q2, W, leptonmass, cospi_sample, 3, params, strucfuncs_sample,strucfuncs_sample);
//         if (ERR == -1){continue;}
//         y[0] = strucfuncs_sample[0][0];
//         y[1] = strucfuncs_sample[1][0];
//         y[2] = strucfuncs_sample[2][0];
       
//         int Nevents=1; //No events now just getting ABC coefficients
//         double events[Nevents];
//         genevents(cospi_sample, y, Nevents, events, A, B, C);

//         //output of events (also total for normalization)
// //        histogram.open(histfnm);
// //        for (int iev = 0; iev < Nevents ; iev++)
// //        {
// //            histogram << events[iev] << "  " << Nevents << endl;
// //        }
// //        histogram.close();

//         for (int j = 0 ; j < cnt ; j++)
//         {
//             double costh = costheta_pi_in[j];
//             double trapfac = 1.;
//             if (j == 0 || j == cnt -1){trapfac = 0.5;}
//             integrated+=0.01*strucfuncs[j][0]*trapfac;
//             chi2 += pow(strucfuncs[j][0] - A*costh*costh - B*costh - C,2)/strucfuncs[j][0];
//         }

//         //output the functions and other stuff
//         distribution.open(funcfnm); 
//         for (int j = 0 ; j < cnt ; j++)
//         {
//             double costh = costheta_pi_in[j];
//             distribution << W << "  " << Q2 << "  " << costh << "  " << strucfuncs[j][0] << "  " << A*costh*costh + B*costh + C <<  "  " << chi2/(cnt-3) <<"  " << integrated << "  " << 2.*A/3. + 2*C << endl;
//         }
//         distribution.close();

//     }
// }
// return 0;
// }
