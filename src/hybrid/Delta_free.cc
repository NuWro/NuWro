#include "Delta_free.h"
#include <math.h>

using namespace std;
#include "Matrix.h"
#include "diracMatrices.h"

#include "Calculations.h"
#include "Constants.h"



// // Oset and Salcedo medium modification of the delta width
void OSMM( int medmod, double s, double r, complex<double> &f_OSMM)
{
    
    f_OSMM = 1.;
    medmod=0; 
    if( medmod == 2 ){
        
        double MDelta2 = MDelta*MDelta;
        
        double W = sqrt(s);
        double EN = ( s - Mpi2 + MN2 )/( 2.*W );

        double qcm2 = EN*EN - MN2;
        double qcm = sqrt(qcm2);
        
        // free decay width
        // // // Leitner et al., Praet's Thesis, ...      
        double Delta_Width = pow( qcm, 3)*( MN + EN )*( fDelta*fDelta )/( 12*Pi * Mpi2*W );
        
        
// // // // Oset and Salcedo recipe:
        
      double q_tilde, I_1, I_2;
	q_tilde = qcm/225;

      // Get I_1 and I_2 from J. Nieves et al., NPA 554, pg. 576 
      if( q_tilde > 1. ){
	I_1 = 1. - 1./( 5*q_tilde*q_tilde );
	I_2 = 1. - 3./( 5*q_tilde*q_tilde ) - 4./( 21*pow(q_tilde,6)) + 18./( 35*pow(q_tilde,4));
      }
      else{
	I_1 = 1. + q_tilde - 1. - 0.2*pow(q_tilde,3);
	I_2 = 1. - 1. + 33./35*q_tilde - 23./105*pow(q_tilde,3);
      }

      // Pauli-blocking corrected width  
      double Pauli_Width = Delta_Width*( I_1 + I_2 )/2;

      double Tpi = ( s - Mpi2 - MN2 )/( 2*MN ) - Mpi;

      double C_QE, C_A2, C_A3;
      double alpha, beta, gamma;
      
      // along Oset and Salcedo, NPA 468, 631 (1987) and J. Nieves (see before)
      if( Tpi < 85. ){
          
        double Tpi_Mpi = Tpi/Mpi;
        double Tpi_Mpi2 = Tpi_Mpi*Tpi_Mpi;

	C_QE = Tpi_Mpi*( 20.2 - 8.58*Tpi_Mpi + 0.702*Tpi_Mpi2 );
	C_A2 = 1.06*Tpi_Mpi2 - 6.64*Tpi_Mpi + 22.66;
	C_A3 = 3.7*Tpi/85;
	alpha = 1 + Tpi_Mpi*( -0.309 - 0.315*Tpi_Mpi + 0.151*Tpi_Mpi2 );
	beta = -0.038*Tpi_Mpi2 + 0.204*Tpi_Mpi + 0.613;
	gamma = 1 + Tpi_Mpi*( 0.984 - 0.512*Tpi_Mpi + 0.1*Tpi_Mpi2 );
        
      }else{
	if( Tpi < 315. ){
            
          double Tpi_Mpi = Tpi/Mpi;
          double Tpi_Mpi2 = Tpi_Mpi*Tpi_Mpi;
            
	  C_QE = -5.19*Tpi_Mpi2 + 15.35*Tpi_Mpi + 2.06;
	  C_A2 = 1.06*Tpi_Mpi2 - 6.64*Tpi_Mpi + 22.66;
	  C_A3 = -13.46*Tpi_Mpi2 + 46.17*Tpi_Mpi - 20.34;
	  alpha = 0.382*Tpi_Mpi2 - 1.322*Tpi_Mpi + 1.466;
	  beta = -0.038*Tpi_Mpi2 + 0.204*Tpi_Mpi + 0.613;
	  gamma = 2*beta;
          
	}else{
            
          double x315_Mpi = 315./Mpi;
          double x315_Mpi2 = x315_Mpi*x315_Mpi;

	  C_QE =  -5.19 *x315_Mpi2 + 15.35*x315_Mpi + 2.06;
	  C_A2 =   1.06 *x315_Mpi2 -  6.64*x315_Mpi + 22.66;
	  C_A3 = -13.46 *x315_Mpi2 + 46.17*x315_Mpi - 20.34;
	  alpha =  0.382*x315_Mpi2 - 1.322*x315_Mpi + 1.466;
	  beta =  -0.038*x315_Mpi2 + 0.204*x315_Mpi + 0.613;
	  gamma =  2*beta;
          
	}
      }

      
// //  density profile:
      double rho;
      double Anum = 12.;
      double rho0 = 0.17; //!fm^-3
      double c = 1.07*pow(Anum,1/3.); //!fm
      double a = 0.54; //!fm
        rho = rho0 / ( 1 + exp((r-c)/a) ); // fm^-3
        
      double rho_rho0 = rho/rho0;
      double QE, A_2, A_3;
	QE = 2.*C_QE*pow(rho_rho0,alpha);
	A_2 = 2.*C_A2*pow(rho_rho0,beta);
	A_3 = 2.*C_A3*pow(rho_rho0,gamma);

        
    // fully corrected decay width, + pion-related part
      double Corr_Width, Corr_Width_Pion;
	Corr_Width = Pauli_Width + QE + A_2 + A_3;
	Corr_Width_Pion = Pauli_Width + QE;

      double Delta_Mass= MDelta;
        // Delta_Mass = MDelta+ 40.*rho_rho0;
        // Delta_Mass = MDelta + 30.;

      // rescale decay coupling to pion part of MM Delta width
      // Notice:    fDelta_scaled = fDelta * sqrt( Corr_Width_Pion / Delta_Width ); //This is the fDelta that should enter in the N-pi_Delta vertex. However, to avoid modify that vertex, we define the quantity "scale = fDelta_scaled / fDelta" and multiply by the PROPAGATOR.
      double scale = sqrt( Corr_Width_Pion / Delta_Width );

      complex<double> Fac_Prop_FF_free = -1./( s - MDelta2 + I*MDelta*Delta_Width );
      
      complex<double> Fac_Prop_FF = scale * (-1.)/ ( s - pow(Delta_Mass,2) + I*Delta_Mass*Corr_Width );

      f_OSMM = Fac_Prop_FF/Fac_Prop_FF_free ;
        
    } //medmod==2
    
}

void Delta_ff(double s, double u, double &Deltaff){

  //Modif, same as in my thesis: cutoff is bigger than in RGJ
  //notice cutoff and VRfac are the same unlike in VR paper, but according to RGJ paper

  double cut_off = 1200.;

  double Lam_piND = cut_off; 
  
  double MDelta2 = MDelta*MDelta;
  double Lam_piND4 = pow( Lam_piND ,4 );
  double VRfac= Lam_piND4;

  
  double Fgauss_u = exp( -pow( u - MDelta2, 2) / Lam_piND4 ); // Gaussian
  double Dipole_u = VRfac / ( pow( u - MDelta2, 2) + VRfac ); // "Dipole"
  double FGaDi_u = Fgauss_u * Dipole_u;

  double Fgauss_s = exp( -pow( s - MDelta2, 2) / Lam_piND4 ); // Gaussian
  double Dipole_s = VRfac / ( pow( s - MDelta2, 2) + VRfac ); // "Dipole"
  double FGaDi_s = Fgauss_s * Dipole_s;

  Deltaff = FGaDi_s + FGaDi_u - FGaDi_s*FGaDi_u;

}

/***************************************************************************
                          THE -W N Delta- VERTEX
***************************************************************************/


void Gamma_WNDelta(int cross, int process,  double Qsq, double s, double Q[], double kResonance[], double ki[], Matrix WNDelta[][4] )
{

double Wrec2 = s;

// // // // PHASES (arXiv:1510.06266v2 [hep-ph] 23 Dec 2015)// // // //
complex<double> PhiV, PhiA;
  PhiV = 1., PhiA = 1.;

int Phases = 1; //Phases == 1 to include phases (and, partially, recover unitarity from LAR paper)
if(Phases == 1 && cross != 1){ 

  double WGeV = sqrt(Wrec2)/1000;    
      
  double Q2GeV = Qsq/1.E6;

  if(Q2GeV>2.5){Q2GeV=2.5;}
  if(WGeV<1.1){WGeV=1.1;}
  if(WGeV>1.4){WGeV=1.4;}
  
  double ow = (WGeV -1.0779);
  
  double ePhiV = 5*ow*( 8.3787 + (2.7315-25.5185*ow)/(0.05308416+pow(0.62862-5*ow,2)) + 301.925*ow - 985.80*pow(ow,2) + 862.025*pow(ow,3) ) * ( pow(1+0.14163*Q2GeV,-2) + (0.066192 + ow*(-0.34057 + 1.631475*ow))*Q2GeV );
      
  PhiV = cos(ePhiV*Pi/180) + I*sin(ePhiV*Pi/180);

  double ePhiA_fitA = 5*ow*(5.2514 + (2.9102-26.5085*ow)/(0.0531901969+pow(0.63033-5*ow,2)) + 266.565*ow - 814.575*pow(ow,2) + 624.05*pow(ow,3) ) * ( pow(1+0.088539*Q2GeV,-2) + (0.026654 + ow*(-1.17305 + 3.66475*ow))*Q2GeV );

  // double ePhiA_fitB = 5*ow*(4.9703 + (2.929-26.6295*ow)/(0.0531256401+pow(0.63051-5*ow,2)) + 264.27*ow - 798.525*pow(ow,2) + 598.85*pow(ow,3) ) * ( pow(1+0.10152*Q2GeV,-2) + (0.041484 + ow*(-1.20715 + 3.7545*ow))*Q2GeV );

  PhiA = cos(ePhiA_fitA*Pi/180.) + I*sin(ePhiA_fitA*Pi/180.);

}
// // // // // // // // // // // // // // // // // // // // // // // //


  /*  
      Now, we get the form factor values, along Lalakulich (PRD 71, 074003) and Leitner (Dipl. Thesis).  
  */
  

  // Use fit by Lalakulich et al. PRD74 014009 (2006)
  
  complex<double> C3_V, C4_V, C5_V;
  C3_V = ( ((2.13)/(pow(1+Qsq/705600,2)))*((1.)/(1+Qsq/(4.*705600.))) ) * PhiV;
  C4_V = -C3_V*1.51/2.13; 
  C5_V = ( ((0.48)/(pow(1+Qsq/705600,2)))*((1.)/(1+Qsq/(0.776*705600.)))) * PhiV;
  
// Amaro 2005 form factors
//   double C3_V, C4_V, C5_V;  
//   C3_V = 2.05/pow((1.+Qsq/540000.),2);
//   C4_V = -C3_V*MN/MDelta;  // we take the proton-mass here, just for simplicity reasons
//   C5_V = 0;
  

  complex<double> C3_A, C4_A, C5_A, C6_A;
  C3_A=0.;
  C4_A=0.;
  C5_A=0.;
  C6_A=0.;

if( process != 0 ){
  
  if(Phases==0){
    C5_A = ((1.2)/(pow(1+Qsq/1102500,2)))*((1.)/(1+Qsq/(3*1102500)));
    C3_A = 0.;
    C4_A = (-1.)/(4)*C5_A;
    C6_A = (MN2)/(Qsq+pow(Mpi,2))*C5_A;  // see remark, following C4_V
  }
// // (arXiv:1510.06266v2 [hep-ph] 23 Dec 2015)  
  else if(Phases==1){ 
    /*FIT A:*/ C5_A = ( (1.12)/(pow(1+Qsq/910116,2)) ) * PhiA ;  
    C3_A = 0.;
    C4_A = (-1.)/(4)*C5_A;
    C6_A = 0;
  }
  
}  
else if( process == 2 ){
    double qw_v = 1.- 2.*Sin2W;
    
    C3_V = qw_v*C3_V;
    C4_V = qw_v*C4_V; 
    C5_V = qw_v*C5_V;
    
}
  
  /*
    In declaring the necessary arrays of Matrices, we distinguish between the vector part and the axial-
    vector part of the -W N Delta- vertex. 
  */

  Matrix NDelta_Vect[4][4];
  Matrix NDelta_Ax[4][4];

  Matrix q_slash;       //ATTENTION this quantity is not the same for delta_cross and delta, we pass -Q[] and Q[], respectively.
    q_slash.M[0][0]=Q[0]       , q_slash.M[0][1]=0.         , q_slash.M[0][2]=-Q[3]        , q_slash.M[0][3]=-Q[1]+I*Q[2],
    q_slash.M[1][0]=0.         , q_slash.M[1][1]=Q[0]       , q_slash.M[1][2]=-Q[1]-I*Q[2], q_slash.M[1][3]=Q[3],
    q_slash.M[2][0]=Q[3]       , q_slash.M[2][1]=Q[1]-I*Q[2], q_slash.M[2][2]=-Q[0]       , q_slash.M[2][3]=0.,
    q_slash.M[3][0]=Q[1]+I*Q[2], q_slash.M[3][1]=-Q[3]      , q_slash.M[3][2]=0.          , q_slash.M[3][3]=-Q[0];

  double q_dot_kD = Q[0]*kResonance[0] - Q[1]*kResonance[1] - Q[2]*kResonance[2] - Q[3]*kResonance[3];

  double q_dot_ki = Q[0]*ki[0] - Q[1]*ki[1] - Q[2]*ki[2] - Q[3]*ki[3];

  
  double sign;
  double block4[4][4];

  complex<double> c3v_mn=C3_V/MN, c4v_mn=C4_V/MN2, c5v_mn=C5_V/MN2;
  
  for(int i=0; i<4; i++){
    if(i == 0){sign = 1.;}else{sign = -1.;}

      block4[i][i] = 1./MN2 *(sign*q_dot_kD - Q[i]*kResonance[i]);
      
      NDelta_Vect[i][i] = (c3v_mn  *(sign*q_slash - Q[i]*Gamma[i])) + (C4_V*block4[i][i] + c5v_mn *(sign*q_dot_ki - Q[i]*ki[i]))*Id;
  }

  
  complex<double> temp, temp3;
  double temp4;
  for(int i=0; i<4; i++){
    temp = c3v_mn  *(-Q[i]);
    temp3 = c5v_mn *(-Q[i]);
    temp4 = -Q[i]/MN2;
    for(int j=0; j<4; j++){
       if( j != i ){
	block4[i][j] = temp4*kResonance[j];
	
	NDelta_Vect[i][j] = temp *Gamma[j] + (C4_V*block4[i][j] + temp3*ki[j])*Id;
       }
    }
  }


    complex<double> c6a_mn;
    if( process==1 || process==2 ){
        // C3_A = 0...
        double gmunu;
        for(int i=0; i<4; i++){
            c6a_mn = C6_A/MN2 * Q[i];
            
            for(int j=0; j<4; j++){
                
                if(i == j){gmunu = 1.;}else{gmunu = 0;}
                if(i == j && i != 0){gmunu = -1.;}
                
                NDelta_Ax[i][j] = /*C3_A * block3[i][i] +*/ ( C4_A * block4[i][j] + C5_A * gmunu +  c6a_mn*Q[j] )*Id;
                
            }
            
        }
    }
    


    if(process == 0){
                
        for(int i=0; i<4; i++){
            for(int j=0; j<4; j++){
                
                WNDelta[i][j].M[0][0] = NDelta_Vect[i][j].M[0][2],
                WNDelta[i][j].M[0][1] = NDelta_Vect[i][j].M[0][3],
                WNDelta[i][j].M[0][2] = NDelta_Vect[i][j].M[0][0],
                WNDelta[i][j].M[0][3] = NDelta_Vect[i][j].M[0][1],

                WNDelta[i][j].M[1][0] = NDelta_Vect[i][j].M[1][2],
                WNDelta[i][j].M[1][1] = NDelta_Vect[i][j].M[1][3],
                WNDelta[i][j].M[1][2] = NDelta_Vect[i][j].M[1][0],
                WNDelta[i][j].M[1][3] = NDelta_Vect[i][j].M[1][1],
                
                WNDelta[i][j].M[2][0] = NDelta_Vect[i][j].M[2][2],
                WNDelta[i][j].M[2][1] = NDelta_Vect[i][j].M[2][3],
                WNDelta[i][j].M[2][2] = NDelta_Vect[i][j].M[2][0],
                WNDelta[i][j].M[2][3] = NDelta_Vect[i][j].M[2][1],
                
                WNDelta[i][j].M[3][0] = NDelta_Vect[i][j].M[3][2],
                WNDelta[i][j].M[3][1] = NDelta_Vect[i][j].M[3][3],
                WNDelta[i][j].M[3][2] = NDelta_Vect[i][j].M[3][0],
                WNDelta[i][j].M[3][3] = NDelta_Vect[i][j].M[3][1];
                
            }
        }   
        
    }else{

        for(int i=0; i<4; i++){
            for(int j=0; j<4; j++){

 		              
                WNDelta[i][j].M[0][0] = NDelta_Vect[i][j].M[0][2] + NDelta_Ax[i][j].M[0][0],
                WNDelta[i][j].M[0][1] = NDelta_Vect[i][j].M[0][3] + NDelta_Ax[i][j].M[0][1],
                WNDelta[i][j].M[0][2] = NDelta_Vect[i][j].M[0][0] + NDelta_Ax[i][j].M[0][2],
                WNDelta[i][j].M[0][3] = NDelta_Vect[i][j].M[0][1] + NDelta_Ax[i][j].M[0][3],

                WNDelta[i][j].M[1][0] = NDelta_Vect[i][j].M[1][2] + NDelta_Ax[i][j].M[1][0],
                WNDelta[i][j].M[1][1] = NDelta_Vect[i][j].M[1][3] + NDelta_Ax[i][j].M[1][1],
                WNDelta[i][j].M[1][2] = NDelta_Vect[i][j].M[1][0] + NDelta_Ax[i][j].M[1][2],
                WNDelta[i][j].M[1][3] = NDelta_Vect[i][j].M[1][1] + NDelta_Ax[i][j].M[1][3],
                
                WNDelta[i][j].M[2][0] = NDelta_Vect[i][j].M[2][2] + NDelta_Ax[i][j].M[2][0],
                WNDelta[i][j].M[2][1] = NDelta_Vect[i][j].M[2][3] + NDelta_Ax[i][j].M[2][1],
                WNDelta[i][j].M[2][2] = NDelta_Vect[i][j].M[2][0] + NDelta_Ax[i][j].M[2][2],
                WNDelta[i][j].M[2][3] = NDelta_Vect[i][j].M[2][1] + NDelta_Ax[i][j].M[2][3],
                
                WNDelta[i][j].M[3][0] = NDelta_Vect[i][j].M[3][2] + NDelta_Ax[i][j].M[3][0],
                WNDelta[i][j].M[3][1] = NDelta_Vect[i][j].M[3][3] + NDelta_Ax[i][j].M[3][1],
                WNDelta[i][j].M[3][2] = NDelta_Vect[i][j].M[3][0] + NDelta_Ax[i][j].M[3][2],
                WNDelta[i][j].M[3][3] = NDelta_Vect[i][j].M[3][1] + NDelta_Ax[i][j].M[3][3];
                
            }
        }
                
    }
            

}


/***********************************************************************
                         THE DELTA PROPAGATOR
***********************************************************************/
void S_Delta( int medmod, int cross, double W2, double kResonance[], Matrix kRSlash, Matrix Delta[][4] )
// void S_Delta( int medmod, int cross, const double &s, const double &u, double kResonance[], Matrix Delta[][4] )
{
  
  double MDelta2 = MDelta*MDelta;

  complex<double> Fac_Prop_FF;
  Matrix Prop_FF;
  double Delta_Width;


  //   // Definition of the energy-dependent Delta decay width \Gamma(W) "Delta_Width"

  if( cross == 1 ){ //the delta-decay width is zero (this is always the case in the Delta-cross pole)

    Fac_Prop_FF = -1./( W2 - MDelta2 );

  }
  if( cross == 0 ){ // we compute the delta-decay width

    double W = sqrt(W2);

    double EN = ( W2 - Mpi2 + MN2 )/( 2.*W );

    double qcm = sqrt(EN*EN - MN2);

    // free decay width
    // // // Leitner et al., Praet's Thesis, ...      
    Delta_Width = pow( qcm, 3)*( MN + EN )*( fDelta*fDelta )/( 12*Pi * Mpi2*W );

    // // // // // // // // // // // //       
    // // // Hernandez et al. 2007,
    //       Delta_Width = ( ( fDelta*fDelta* MN * pow( qcm, 3) )/( 6*Pi * Mpi2*W ) );
    // // // // // // // // // // // // 
    
    Fac_Prop_FF = -1./( W2 - MDelta2 + I*MDelta*Delta_Width );
      
  } // //if cross==0

  Prop_FF = (Fac_Prop_FF/3.) *( kRSlash + MDelta*Id ); //The factor 1/3. here is from the Rarita propagator

// // // // // // // // // // // // // // // // // // // // // // // // // // // // // //   
  Matrix Rarita[4][4];

  double a0=2. - 2./MDelta2 *pow(kResonance[0],2), a1=-2 - 2./MDelta2 *pow(kResonance[1],2), a2 = -2. - 2./MDelta2 *pow(kResonance[2],2), a3 = -2. - 2./MDelta2 *pow(kResonance[3],2);

Matrix Sym, Asym;
  for(int i=0; i<4; i++){
  
    Rarita[0][0].M[i][i] = a0 ; 
    Rarita[1][1].M[i][i] = a1 ; 
    Rarita[2][2].M[i][i] = a2 ; 
    Rarita[3][3].M[i][i] = a3 ; 
      
    for(int j=i+1; j<4; j++){
    
    Sym = (-2./MDelta2*kResonance[i]*kResonance[j])*Id;
    Asym= ( (kResonance[i]/MDelta)*Gamma[j] - (kResonance[j]/MDelta)*Gamma[i]);
	Rarita[i][j] = mGamma_munu[i][j] + Sym + Asym;
    Rarita[j][i] = mGamma_munu[j][i] + Sym - Asym;
        

    }
  }

  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){

      Delta[i][j] = Prop_FF*( Rarita[i][j] );

    }
  }

}


/***************************************************************************
                          THE -Delta pi N- VERTEX
***************************************************************************/

/***************************************************************************
The vertex used by J. Weda (PhD thesis) is taken. 
Only used if Pascalutsa = 1;
***************************************************************************/

void Gamma_DeltaNpi( double s, double fDeltaNpi, double kResonance[], double kpi[], Matrix DeltaNpi[] ){

  int epsilon[4][4][4][4];
     for(int i=0; i<4; i++){
     for(int j=0; j<4; j++){
     for(int k=0; k<4; k++){
     for(int l=0; l<4; l++){
       epsilon[i][j][k][l] = 0;
     }
     }
     }
     }
      epsilon[0][1][2][3] = -1; epsilon[0][1][3][2] = 1;  
      epsilon[0][2][1][3] = 1; epsilon[0][2][3][1] = -1; 
      epsilon[0][3][2][1] = 1; epsilon[0][3][1][2] = -1; 
    
      epsilon[1][0][2][3] = 1; epsilon[1][0][3][2] = -1; 
      epsilon[1][2][0][3] = -1; epsilon[1][2][3][0] = 1; 
      epsilon[1][3][2][0] = -1; epsilon[1][3][0][2] = 1;
      
      epsilon[2][1][0][3] = 1; epsilon[2][1][3][0] = -1;  
      epsilon[2][0][1][3] = -1; epsilon[2][0][3][1] = 1; 
      epsilon[2][3][0][1] = -1; epsilon[2][3][1][0] = 1; 
      
      epsilon[3][0][2][1] = -1; epsilon[3][0][1][2] = 1; 
      epsilon[3][2][0][1] = 1; epsilon[3][2][1][0] = -1; 
      epsilon[3][1][2][0] = 1; epsilon[3][1][0][2] = -1;
      

// // //  covariant vertex       
      
    for(int i=0; i<4; i++){
      for(int j=0; j<4; j++){
	if(j != i){
	  for(int k=0; k<4; k++){
	    if(k != i && k != j){
	      for(int l=0; l<4; l++){
		if( l != i && l != j && l != k ){
		  
		  DeltaNpi[i] = DeltaNpi[i] + (fDeltaNpi)/(Mpi*MDelta)* epsilon[i][j][k][l] *kpi[j]* kResonance[l] * Gamma[k]*Gamma5 ;
		
		}
	      }
	    }
	  }
	}
      }
    }

}

void DP_current( int medmod, int Pascalutsa, int process, int nucleon, int decay, int Helicity, int cross, double Qsq, double s, double u, double Q[], double kResonance[], double ki[], double kpi[], Matrix kRSlash, Matrix Op_delta[] )
{
    double W2 = s;
// // // // NOTATION // // // // 
// process = 1   CC interaction
//   HELICITY == -1 --> W^+ induced 1-pion production 
// 	nucleon = 1 --> proton initial state
// 		decay = 1 or 2 --> p + pi^+
// 	nucleon = 2 --> neutron initial state
// 		decay = 1 --> p + pi^0
// 		decay = 2 --> n + pi^+
//   
//   HELICITY == 1 --> W^- induced 1-pion production 
// 	nucleon = 2 --> neutron initial state
// 		decay = 1 or 2 --> n + pi^-
// 	nucleon = 1 --> proton initial state
// 		decay = 1 --> n + pi^0
// 		decay = 2 --> p + pi^-
// //    
// process = 0   EM interaction --> photon induced 1-pion production 
// 	nucleon = 1 --> proton initial state 
// 		decay = 1 --> p + pi^0
// 		decay = 2 --> n + pi^+
// 	nucleon = 2 --> neutron initial state
// 		decay = 1 --> n + pi^0
// 		decay = 2 --> p + pi^-
// // // // // // // // // // //
  
// // //  Delta 
      double icD;
if( process == 1 ){
  
  if(Helicity == -1){ //neutrino
      if( nucleon == 1 ){ 
	icD = sqrt(3./2.);
      }
      else{ 
	if( decay == 1 ){ 
	  icD = -sqrt(1./3.); 
	}else{ 
	  icD = sqrt(1./6.); 
	} 
      }
  }
  else if(Helicity == 1){ //antineutrino
      if( nucleon == 2 ){ 
	icD = sqrt(3./2.);
      }
      else{ 
	if( decay == 1 ){ icD = sqrt(1./3.); }else{ icD = sqrt(1./6.); } 
      }
  }
  
}
else if( process == 0 || process == 2 ){
  
  if( nucleon == 1 ){
    if( decay == 1 ){ icD = sqrt(1./3.);}else{ icD = -sqrt(1./6.) ;}
  }
  else if( nucleon == 2 ){
    if( decay == 1 ){ icD = sqrt(1./3.);}else{ icD = sqrt(1./6.) ;}
  }
  
}



if( icD != 0 ){
//       // -W N Delta- vertex 
      Matrix WNDelta[4][4];
      Gamma_WNDelta( cross, process, Qsq, s, Q, kResonance, ki, WNDelta );
//       // Delta propagator
      Matrix Delta[4][4];
      S_Delta( medmod, cross, W2, kResonance, kRSlash, Delta );
// //  Delta[i][j] has been defined in CONTRAVARIANT (upper indices) notation but we need it in covariant (down indices)
      for( int i=1; i<4; i++ ){
	  Delta[0][i] = (-1)*Delta[0][i];
	  Delta[i][0] = (-1)*Delta[i][0];
	}
//       // -Delta pi N- vertex   
      Matrix DeltaNpi[4];
      
      if(Pascalutsa == 0){
	for( int i=0; i<4; i++ ){ 
	  DeltaNpi[i] =  (icD* fDelta*sqrt(2.)/Mpi *kpi[i])*Id; 
	}     
      }
      if(Pascalutsa == 1){
      double fDeltaNpi = icD*fDelta*sqrt(2.);
      Gamma_DeltaNpi( W2, fDeltaNpi, kResonance, kpi, DeltaNpi );
      }

    Matrix temp;
    for( int k=0; k<4; k++ ){
        for( int l=0; l<4; l++ ){
            temp = DeltaNpi[k]*Delta[k][l];
            for( int i=0; i<4; i++ ){
        
                Op_delta[i] = Op_delta[i] + temp*WNDelta[l][i];
            
	  }
	}
      }
    
      
  } //icD different than cero
}


void DP_current_pre( int medmod, int Pascalutsa, int process, int nucleon, int decay, int Helicity, int cross, double Qsq, double s, double u, double Q[], double kResonance[], double ki[], double kpi[], Matrix kRSlash,
 const Matrix WNDelta[][4],
 const Matrix PropDelta[][4],
Matrix Op_delta[] )
{
    double W2 = s;
// // // // NOTATION // // // // 
// process = 1   CC interaction
//   HELICITY == -1 --> W^+ induced 1-pion production 
// 	nucleon = 1 --> proton initial state
// 		decay = 1 or 2 --> p + pi^+
// 	nucleon = 2 --> neutron initial state
// 		decay = 1 --> p + pi^0
// 		decay = 2 --> n + pi^+
//   
//   HELICITY == 1 --> W^- induced 1-pion production 
// 	nucleon = 2 --> neutron initial state
// 		decay = 1 or 2 --> n + pi^-
// 	nucleon = 1 --> proton initial state
// 		decay = 1 --> n + pi^0
// 		decay = 2 --> p + pi^-
// //    
// process = 0   EM interaction --> photon induced 1-pion production 
// 	nucleon = 1 --> proton initial state 
// 		decay = 1 --> p + pi^0
// 		decay = 2 --> n + pi^+
// 	nucleon = 2 --> neutron initial state
// 		decay = 1 --> n + pi^0
// 		decay = 2 --> p + pi^-
// // // // // // // // // // //
  
// // //  Delta 
      double icD;
if( process == 1 ){
  
  if(Helicity == -1){ //neutrino
      if( nucleon == 1 ){ 
	icD = sqrt(3./2.);
      }
      else{ 
	if( decay == 1 ){ 
	  icD = -sqrt(1./3.); 
	}else{ 
	  icD = sqrt(1./6.); 
	} 
      }
  }
  else if(Helicity == 1){ //antineutrino
      if( nucleon == 2 ){ 
	icD = sqrt(3./2.);
      }
      else{ 
	if( decay == 1 ){ icD = sqrt(1./3.); }else{ icD = sqrt(1./6.); } 
      }
  }
  
}
else if( process == 0 || process == 2 ){
  
  if( nucleon == 1 ){
    if( decay == 1 ){ icD = sqrt(1./3.);}else{ icD = -sqrt(1./6.) ;}
  }
  else if( nucleon == 2 ){
    if( decay == 1 ){ icD = sqrt(1./3.);}else{ icD = sqrt(1./6.) ;}
  }
  
}


 if( icD == 0 ){
   for( int i=0; i<4; i++ ){
     Op_delta[i] = 0*Id;
   }  
 }

if( icD != 0 ){
//       // -Delta pi N- vertex   
      Matrix DeltaNpi[4];
      
      if(Pascalutsa == 0){
	for( int i=0; i<4; i++ ){ 
	  DeltaNpi[i] =  (icD* fDelta*sqrt(2.)/Mpi *kpi[i])*Id; 
	}     
      }
      if(Pascalutsa == 1){
      double fDeltaNpi = icD*fDelta*sqrt(2.);
      Gamma_DeltaNpi( W2, fDeltaNpi, kResonance, kpi, DeltaNpi );
      }

    Matrix temp;
    for( int k=0; k<4; k++ ){
        for( int l=0; l<4; l++ ){
            temp = DeltaNpi[k]*PropDelta[k][l];
            for( int i=0; i<4; i++ ){
        
                Op_delta[i] = Op_delta[i] + temp*WNDelta[l][i];
            
	  }
	}
      }
    
      
  } //icD different than cero
}


void CDP_current( int medmod, int Pascalutsa, int process, int nucleon, int decay, int Helicity, int cross, double Qsq, double s, double u, double minusQ[], double kResonance[], double kN[], double kpi[], Matrix kRSlash, Matrix Op_delta_cross[] )
{      

    double W2 = u;
    
// // Delta crossed contribution
      double icD;
if( process == 1 ){
  
  if(Helicity == -1){ //neutrino
    
      if( nucleon == 1 ){ 
	icD = sqrt(1./6.);
      }
      else{ 
	if( decay == 1 ){ 
	  icD = sqrt(1./3.); 
	}else{ 
	  icD = sqrt(3./2.); 
	} 
      }
      
  }
  else if(Helicity == 1){ //antineutrino
      if( nucleon == 2 ){ 
	icD = sqrt(1./6.);
      }
      else{ 
	if( decay == 1 ){ 
	  icD = -sqrt(1./3.); 
	}else{ 
	  icD = sqrt(3./2.); 
	} 
      }
  }
  
}
else if( process == 0 || process == 2 ){
  if( nucleon == 1 ){
    if( decay == 1 ){ icD = sqrt(1./3.);}else{ icD = sqrt(1./6.) ;}
  }
  else if( nucleon == 2 ){
    if( decay == 1 ){ icD = sqrt(1./3.);}else{ icD = -sqrt(1./6.) ;}
  }
}




if(icD != 0){
//       // -W N Delta- vertex

      Matrix WNDelta_cross[4][4];
      Gamma_WNDelta( cross, process, Qsq, s, minusQ, kResonance, kN, WNDelta_cross );
//       // Delta propagator
      Matrix Delta_cross[4][4];     
      S_Delta( medmod, cross, W2, kResonance, kRSlash, Delta_cross );
// Delta_cross[i][j] has been defined in CONTRAVARIANT (upper indices) notation but we need it in covariant (down indices)
      for( int i=1; i<4; i++ ){
	  Delta_cross[0][i] = (-1.)*Delta_cross[0][i];
	  Delta_cross[i][0] = (-1.)*Delta_cross[i][0];
	}
// // // // //
//       // -Delta pi N- vertex 
      Matrix DeltaNpi_cross[4];
      if(Pascalutsa == 0){
	for( int i=0; i<4; i++ ){ 
	  DeltaNpi_cross[i] =  ( icD*fDelta*sqrt(2.)/Mpi *kpi[i])*Id; 
	}    
    }
   if(Pascalutsa == 1){
        double fDeltaNpi = icD*fDelta*sqrt(2.);
        Gamma_DeltaNpi( W2, fDeltaNpi, kResonance, kpi, DeltaNpi_cross );
      }	
      
      Matrix block_delta_c[4][4];
      for( int i=0; i<4; i++ ){
	for( int j=0; j<4; j++ ){
	  block_delta_c[i][j] = ConjMatrix( WNDelta_cross[j][i] );
          
                block_delta_c[i][j].M[0][2] = -block_delta_c[i][j].M[0][2];
                block_delta_c[i][j].M[0][3] = -block_delta_c[i][j].M[0][3];
                block_delta_c[i][j].M[1][2] = -block_delta_c[i][j].M[1][2];
                block_delta_c[i][j].M[1][3] = -block_delta_c[i][j].M[1][3];
                
                block_delta_c[i][j].M[2][0] = -block_delta_c[i][j].M[2][0];
                block_delta_c[i][j].M[2][1] = -block_delta_c[i][j].M[2][1];
                block_delta_c[i][j].M[3][0] = -block_delta_c[i][j].M[3][0];
                block_delta_c[i][j].M[3][1] = -block_delta_c[i][j].M[3][1];
                
	}
      }
      

    Matrix temp;
    for( int k=0; k<4; k++ ){
        for( int l=0; l<4; l++ ){
            temp = Delta_cross[l][k]*DeltaNpi_cross[k];
            for( int i=0; i<4; i++ ){
          
                Op_delta_cross[i] = Op_delta_cross[i] + block_delta_c[i][l]*temp;
            
	  }
	}
      }  


    } //if(icD != 0)
    
}

