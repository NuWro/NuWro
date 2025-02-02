#include "NON-Resonant.h"
// /***************************************************************************
//                           THE -W N N- VERTEX
// ***************************************************************************/

void Gamma_WNN( int nucleon, int process, int decay, int cross, double Qsq, double Q[], const Matrix& QSlash, Matrix WNN[] ){
 
  double QsqGeV, xmu, tau;
  QsqGeV = Qsq/1.E6;
  xmu = 2.*MN;
  tau = Qsq/(xmu*xmu);
      
  double DipV, DipA, MA2, MV2;
  MA2 = 1.1025; //MA=1.05 GeV
  MV2 = 0.710649; //MV=0.843 GeV
  DipV = pow( 1+QsqGeV/MV2, 2 ); 
  DipA = pow( 1+QsqGeV/MA2, 2 );  
    
  
    double F1p, F2p, F1n, F2n, GA0, GP0, F1V, F2V;   

    double GEp, GEn, GMp, GMn;
    
// //  GALSTER
    GEp = 1./DipV;
    GMp = 2.792847*GEp;
    GEn = 1.913043*tau*GEp/(1.+5.6*tau);
    GMn = -1.913043*GEp;
// // // // // // 
    
// //  KELLY
// 	double rmup=2.79285, rmun=-1.91304;
// 	
// 	double As=1.70, Bs=3.30;
// 	
// 	double a0=1.0;
// 
// 	double a1Gep=-0.24, b1Gep=10.98, b2Gep=12.82, b3Gep=21.97;
// 
// 	double a1Gmp=0.12, b1Gmp=10.97, b2Gmp=18.86, b3Gmp=6.55;
// 	
// 	double a1Gmn=2.33, b1Gmn=14.72, b2Gmn=24.20, b3Gmn=84.1;
// 
// 	double GD=1./pow((1.+QsqGeV/0.71),2);
// 
// 	GEn=As*tau/(1.+Bs*tau)*GD;
// 	GEp=(a0 + a1Gep*tau)/(1. + b1Gep*tau + b2Gep*pow(tau,2) + b3Gep*pow(tau,3));
// 	GMp=rmup*(a0 + a1Gmp*tau)/(1. + b1Gmp*tau + b2Gmp*pow(tau,2) + b3Gmp*pow(tau,3));
// 	GMn=rmun*(a0 + a1Gmn*tau)/(1. + b1Gmn*tau + b2Gmn*pow(tau,2) + b3Gmn*pow(tau,3));
// // // // // //     
    
    F1p = (GEp+tau*GMp)/(1+tau);
    F2p = (GMp-GEp)/(1+tau);
    F1n = (GEn+tau*GMn)/(1+tau);
    F2n = (GMn-GEn)/(1+tau);
  
    GA0 = gA/ DipA ;

  double F1, F2; //Set to correct values in the following
  F1=0.;
  F2=0.;
  
  if( process == 0 ){
    if( nucleon == 1 && decay == 1){      
      F1 = F1p;
      F2 = F2p;
    }
    if( nucleon == 1 && decay == 2){      
      if( cross == 0 ){
	F1 = F1p;
	F2 = F2p;
      }
      if( cross == 1 ){
	F1 = F1n;
	F2 = F2n;
      }
    }    
    if( nucleon == 2 && decay == 1){
      F1 = F1n;
      F2 = F2n;
    }
    if( nucleon == 2 && decay == 2){
      if( cross == 0){
	F1 = F1n;
	F2 = F2n;
      }
      if( cross == 1 ){
	F1 = F1p;
	F2 = F2p;
      }
    }
  }
  
if( process == 2 ){
    
  double F1s = 0., F2s = 0.;
  double GAs = 0.;
  
  double wF1p = 0.5*(QWeak*F1p - F1n - F1s); 
  double wF2p = 0.5*(QWeak*F2p - F2n - F2s); 

  double wF1n = 0.5*(QWeak*F1n - F1p - F1s); 
  double wF2n = 0.5*(QWeak*F2n - F2p - F2s); 

  double wGAp = 0.5*(GA0 - GAs);
  double wGAn = 0.5*(-GA0 - GAs);
  
  if( nucleon == 1 && decay == 1){      
    F1 = wF1p;
    F2 = wF2p;
    
    GA0 = wGAp;
  }
  if( nucleon == 1 && decay == 2){      
    if( cross == 0 ){
      F1 = wF1p;
      F2 = wF2p;
    
      GA0 = wGAp;
    }
    if( cross == 1 ){
      F1 = wF1n;
      F2 = wF2n;
    
      GA0 = wGAn;
    }
  }    
  if( nucleon == 2 && decay == 1){
    F1 = wF1n;
    F2 = wF2n;
    
    GA0 = wGAn;
  }
  if( nucleon == 2 && decay == 2){
    if( cross == 0){
      F1 = wF1n;
      F2 = wF2n;
    
      GA0 = wGAn;
    }
    if( cross == 1 ){
      F1 = wF1p;
      F2 = wF2p;
    
      GA0 = wGAp;
    }
  }
}
// // // //   

// // 
  if( process == 1 ){
    F1 = (F1p - F1n);
    F2 = (F2p - F2n);    
  }
// // // 

  GP0 = GA0*pow(MN,2)/(Qsq + pow(Mpi,2)); 
    
    
  Matrix WNN_V[4];
  Matrix WNN_A[4];

// // // // // // Vector part // // // // // // // // // // // // // // 
  for(int i=0; i<4; i++){
    WNN_V[i] = F1*Gamma[i] - F2/(xmu)*0.5*(Gamma[i]*QSlash-QSlash*Gamma[i]);
  }
 
// // // // // // Axial part // // // // // // // // // // // // // // 
  
  if( process == 1 || process == 2 ){
    Matrix QSlash_5;
    QSlash_5 = QSlash*Gamma5;
    for(int i=0; i<4; i++){
      WNN_A[i] = GA0*Gamma_mu5[i] + (GP0/pow(MN,2)*Q[i])*QSlash_5;
    }
  }

//       


  for(int i=0; i<4; i++){

    if( process == 0 ){
      WNN[i] = WNN_V[i] ;
    }
    if( process == 1 || process == 2 ){    
      WNN[i] = WNN_V[i] - WNN_A[i] ;
    }
    
  }

}


/***********************************************************************
                         THE Nucleon PROPAGATOR
***********************************************************************/


void S_Nprop( double W2, const Matrix& kSlash, Matrix &Nprop ){
  
  Nprop = 1. /( W2 - MN2 )*( kSlash + MN*Id );
      
}


void NP_current( int process, int nucleon, int decay, int Helicity, int cross, double Qsq, double W2, double Q[], const Matrix& QSlash, const Matrix& sSlash,const Matrix& kpiSlash, Matrix Op_NP[] )
{
  
// // // // // // // // ISOSPIN FACTORS // // // // // // // 
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
// // // // // // // // // // // // // // // // // // // //
      double icNP=0; 
      
if( cross == 0 ){
  if( process == 0 || process == 2 ){
    if( nucleon==1 ){ 
      if( decay==1 ){ icNP = sqrt(1./2.); } 
      if( decay==2 ){ icNP = 1.; } 
    }
    if( nucleon == 2 ){ 
      if( decay==1 ){ icNP = -sqrt(1./2.); } 
      if( decay==2 ){ icNP = 1.; } 
     }
  }
// // // 
  if( process == 1 ){
    if( Helicity == -1 ){
      if( nucleon == 1 ){ icNP = 0.; } 
      if( nucleon == 2 ){ 
	if( decay == 1 ){ icNP = sqrt(1./2.); } 
	if( decay == 2 ){ icNP = 1.; } 
      }
   }
    if( Helicity == 1 ){
      if( nucleon == 2 ){ icNP = 0.; } 
      if( nucleon == 1 ){ 
	if( decay == 1 ){ icNP = -sqrt(1./2.); } 
	if( decay == 2 ){ icNP = 1.; } 
      }
   }
  }
}
// // // 
// // // 
if( cross == 1 ){
  if( process == 0 || process == 2 ){
    if( nucleon == 1 ){ 
      if( decay == 1 ){ icNP = sqrt(1./2.); } 
      if( decay == 2 ){ icNP = 1.; } 
    }
    if( nucleon == 2 ){ 
      if( decay == 1 ){ icNP = -sqrt(1./2.); } 
      if( decay == 2 ){ icNP = 1.; } 
     }
  }
// // // 
  if( process == 1 ){
    if( Helicity == -1 ){
      if( nucleon == 1 ){ icNP = 1.; } 
      if( nucleon == 2 ){ 
	if( decay == 1 ){ icNP = -sqrt(1./2.); } 
	if( decay == 2 ){ icNP = 0.; } 
      }
    }
    if( Helicity == 1 ){
      if( nucleon == 2 ){ icNP = 1.; } 
      if( nucleon == 1 ){ 
	if( decay == 1 ){ icNP = sqrt(1./2.); } 
	if( decay == 2 ){ icNP = 0.; } 
      }
    }
  }
}

if( icNP == 0 ){
  for( int i=0; i<4; i++ ){
    Op_NP[i] = 0.*Id;
  }
}
else{
//       // -W N N- vertex      
      Matrix WNN[4];
      Gamma_WNN( nucleon, process, decay, cross, Qsq, Q,QSlash, WNN ); 
//       // NUCLEON PROPAGATOR (NP)
      Matrix Nprop;
      S_Nprop( W2, sSlash, Nprop );
//       // -N pi N- vertex 
      Matrix NNpi;
      double fact = icNP * (-gA)/(fpi*sqrt(2));
      NNpi = fact*(kpiSlash*Gamma5); 
      
      Matrix block_NP;
      block_NP = NNpi*Nprop;

      for( int i=0; i<4; i++ ){
	Op_NP[i] = block_NP*WNN[i];
      } 
  }
}


void CNP_current( int process, int nucleon, int decay, int Helicity, int cross, double Qsq, double W_cross2, double Q[], const Matrix& QSlash, const Matrix& uSlash, const Matrix& kpiSlash, Matrix Op_NP_cross[] )
{      

// // // // // // // // ISOSPIN FACTORS // // // // // // // 
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
// 		decay = 2 --> p + pi^+
// //    
// process = 0   EM interaction --> photon induced 1-pion production 
// 	nucleon = 1 --> proton initial state 
// 		decay = 1 --> p + pi^0
// 		decay = 2 --> n + pi^+
// 	nucleon = 2 --> neutron initial state
// 		decay = 1 --> n + pi^0
// 		decay = 2 --> p + pi^-
// // // // // // // // // // //
// // // // // // // // // // // // // // // // // // // // 
  
      double icNP=0; 
if( cross == 0 ){
  if( process == 0 || process == 2 ){
    if( nucleon == 1 ){ 
      if( decay == 1 ){ icNP = sqrt(1./2.); } 
      if( decay == 2 ){ icNP = 1.; } 
    }
    if( nucleon == 2 ){ 
      if( decay == 1 ){ icNP = -sqrt(1./2.); } 
      if( decay == 2 ){ icNP = 1.; } 
     }
  }
// // // 
  if( process == 1 ){
    if( Helicity == -1 ){
      if( nucleon == 1 ){ icNP = 0.; } 
      if( nucleon == 2 ){ 
	if( decay == 1 ){ icNP = sqrt(1./2.); } 
	if( decay == 2 ){ icNP = 1.; } 
      }
   }
    if( Helicity == 1 ){
      if( nucleon == 2 ){ icNP = 0.; } 
      if( nucleon == 1 ){ 
	if( decay == 1 ){ icNP = -sqrt(1./2.); } 
	if( decay == 2 ){ icNP = 1.; } 
      }
   }
  }
}
// // // 
// // // 
if( cross == 1 ){
  if( process == 0 || process == 2 ){
    if( nucleon == 1 ){ 
      if( decay == 1 ){ icNP = sqrt(1./2.); } 
      if( decay == 2 ){ icNP = 1.; } 
    }
    if( nucleon == 2 ){ 
      if( decay == 1 ){ icNP = -sqrt(1./2.); } 
      if( decay == 2 ){ icNP = 1.; } 
     }
  }
// // // 
  if( process == 1 ){
    if( Helicity == -1 ){
      if( nucleon == 1 ){ icNP = 1.; } 
      if( nucleon == 2 ){ 
	if( decay == 1 ){ icNP = -sqrt(1./2.); } 
	if( decay == 2 ){ icNP = 0.; } 
      }
    }
    if( Helicity == 1 ){
      if( nucleon == 2 ){ icNP = 1.; } 
      if( nucleon == 1 ){ 
	if( decay == 1 ){ icNP = sqrt(1./2.); } 
	if( decay == 2 ){ icNP = 0.; } 
      }
    }
  }
}


if( icNP == 0 ){
  for( int i=0; i<4; i++ ){
    Op_NP_cross[i] = 0.*Id;
  } 
}
else{
//       // -W N N- vertex      
      Matrix WNN_cross[4];
      Gamma_WNN( nucleon, process, decay, cross, Qsq, Q, QSlash, WNN_cross ); 
//       // NUCLEON PROPAGATOR (NP)
      Matrix Nprop_cross;
      S_Nprop( W_cross2, uSlash, Nprop_cross );
//       // -N pi N- vertex 
      Matrix NNpi_cross;  
      double fact = icNP * (-gA)/(fpi*sqrt(2.));
      NNpi_cross = fact * (kpiSlash*Gamma5); 
      
      Matrix block_NP_c;
      block_NP_c = Nprop_cross*NNpi_cross;
      
      for( int i=0; i<4; i++ ){
	Op_NP_cross[i] = WNN_cross[i]*block_NP_c;
      }
  }
  
}


