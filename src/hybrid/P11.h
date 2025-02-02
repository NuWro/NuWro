#ifndef P11_H
#define P11_H


void P11_ff( double s, double u, double &P11ff)
{
    
  double cut_off = 1200.;
  double Lam_piND = cut_off; 
  
  double MP112=pow(MP11,2);
  double Lam_piND4 = pow( Lam_piND ,4 );
  
  double Fgauss_u = exp( -pow( u - MP112, 2) / Lam_piND4 ); // Gaussian
  double Dipole_u = Lam_piND4 / ( pow( u - MP112, 2) + Lam_piND4 ); // "Dipole"
  double FGaDi_u = Fgauss_u * Dipole_u;

  double Fgauss_s = exp( -pow( s - MP112, 2) / Lam_piND4 ); // Gaussian
  double Dipole_s = Lam_piND4 / ( pow( s - MP112, 2) + Lam_piND4 ); // "Dipole"
  double FGaDi_s = Fgauss_s * Dipole_s;

  P11ff = FGaDi_s + FGaDi_u - FGaDi_s*FGaDi_u;

}

// /***************************************************************************
//                           THE -W N R- VERTEX
// ***************************************************************************/

void Gamma_WNP11( int nucleon, int process, int decay, int cross, double Qsq, double Q[], double kResonance[], double s, Matrix QSlash, Matrix WNR[] ){
 
  
  double QsqGeV, xmu, tau;
  QsqGeV = Qsq/1.E6;
  xmu = MN + MP11;
    
  double DipV, DipA, MA2, MV2;
  MA2 = 1.1025; //MA=1.05 GeV
  MV2 = 0.7100; //MV=0.84 GeV
  DipV = pow( 1+QsqGeV/MV2, 2 ); 
  DipA = pow( 1+QsqGeV/MA2, 2 );  
    
  
  double F1p, F2p, F1n, F2n, GA0, GP0, F1V, F2V;  

	
   // axial coupling from PCAC, Q2 dependence as in Lalakulich 06
    GA0 = 0.51 / ( DipA * (1. + QsqGeV/(3.*MA2)) ) ;
 
// Vector form factors
// Hernandez result, with W fixed to resonance mass 
// Agrees well with MAID07 for isovector current, but slightly different for isoscalar!
    double W = MP11;
    F1p = -5.7/(DipV*(1. + QsqGeV/(1.4*MV2)));
    F2p = -0.64/DipV*( 1. - 2.47*log(1.+QsqGeV/1.0) );
    
    F1V = ( F1p*(pow(MN+W,2)+ 5./3.*Qsq) + 2./3*F2p*(MN+W)*xmu  ) / ( pow(MN + W,2) + Qsq  );
    F2V = ( F2p*(5.*pow(MN+W,2)+ 3.*Qsq)*xmu + 2.*F1p*Qsq*(MN+W)  ) / ( 3.*xmu*(pow(MN + W,2) + Qsq)  );
    
    F1n = F1p - F1V;
    F2n = F2p - F2V;
 
//////////////////////////////////////////////////
/* OTHER CHOICES FOR Vector FF  
  int Vff = 1;
  
  if( Vff == 1 ){
// // // Lalakulich 2006   // oppositer signs wrt Lalakulich2006
    F1p = -2.3/( DipV * (1.+ QsqGeV/(4.3*MV2)) );
    F2p = -0.76/DipV*( 1.- 0.28*log(1.+QsqGeV/1.0) );
    F1n = -F1p;
    F2n = -F2p;
  
    GA0 = 0.51 / ( DipA * (1. + QsqGeV/(3.*MA2)) ) ;
// // // // // // // // // // // //  
    
  }
  else if( Vff == 2 ){
   double W = sqrt(s);
// // // Hernandez 2008  WARNING : F1V, F2V and therefore F1n, F2n depend on W. However, W is only well define for diract terms (s-channel), for cross terms (u-channels) W2 may be negative and W would be imaginary!! Therefore: don't use this prescription in cross!!
// 
    F1p = -5.7/(DipV*(1. + QsqGeV/(1.4*MV2)));
    F2p = -0.64/DipV*( 1. - 2.47*log(1.+QsqGeV/1.0) );
    
    F1V = ( F1p*(pow(MN+W,2)+ 5./3.*Qsq) + 2./3*F2p*(MN+W)*xmu  ) / ( pow(MN + W,2) + Qsq  );
    F2V = ( F2p*(5.*pow(MN+W,2)+ 3.*Qsq)*xmu + 2.*F1p*Qsq*(MN+W)  ) / ( 3.*xmu*(pow(MN + W,2) + Qsq)  );
    
    F1n = F1p - F1V;
    F2n = F2p - F2V;
  
    GA0 = 0.63 / DipA;
// // // // // // // // // // // //   
  }
/////////////////////////////////////
END OTHER VFF
////////////////////////////////////////
/// */
  
  
double F1, F2;

if( process == 0 ){
  if( nucleon == 1 && decay == 1){      
    F1 = F1p;
    F2 = F2p;
  }
  else if( nucleon == 1 && decay == 2){      
    if( cross == 0 ){
      F1 = F1p;
      F2 = F2p;
    }
    else if( cross == 1 ){
      F1 = F1n;
      F2 = F2n;
    }
  }    
  else if( nucleon == 2 && decay == 1){
    F1 = F1n;
    F2 = F2n;
  }
  else if( nucleon == 2 && decay == 2){
    if( cross == 0){
      F1 = F1n;
      F2 = F2n;
    }
    else if( cross == 1 ){
      F1 = F1p;
      F2 = F2p;
    }
  }
}
// // // // 
else if( process == 2 ){
    
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
  else if( nucleon == 1 && decay == 2){      
    if( cross == 0 ){
      F1 = wF1p;
      F2 = wF2p;
    
      GA0 = wGAp;
    }
    else if( cross == 1 ){
      F1 = wF1n;
      F2 = wF2n;
    
      GA0 = wGAn;
    }
  }    
  else if( nucleon == 2 && decay == 1){
    F1 = wF1n;
    F2 = wF2n;
    
    GA0 = wGAn;
  }
  else if( nucleon == 2 && decay == 2){
    if( cross == 0){
      F1 = wF1n;
      F2 = wF2n;
    
      GA0 = wGAn;
    }
    else if( cross == 1 ){
      F1 = wF1p;
      F2 = wF2p;
    
      GA0 = wGAp;
    }
  }
}
// // // // // // // //   

// // // // //   
if( process == 1 ){
    F1 = (F1p - F1n);
    F2 = (F2p - F2n);
}
// // // // //  

GP0 = xmu/(Qsq + pow(Mpi,2)) * GA0; //MeV^-1
    
  
  
if( process == 1 || process == 2 ){
    for(int i=0; i<4; i++){
        WNR[i] = ( F1/(pow(xmu,2))*(Qsq*Gamma[i] + Q[i]*QSlash) - 0.5*F2/(xmu)*(Gamma[i]*QSlash-QSlash*Gamma[i]) ) + ( GA0*Gamma_mu5[i] + (GP0* Q[i])*Gamma5 ) ;
    }
}
// // // // // // // // // // // // // // // // // // // // // // // // // 
else{
    for(int i=0; i<4; i++){

      WNR[i] = F1/(pow(xmu,2))*(Qsq*Gamma[i] + Q[i]*QSlash) - 0.5*F2/(xmu)*(Gamma[i]*QSlash-QSlash*Gamma[i]) ;
      
    }
  }

}


/***********************************************************************
                         THE P11 PROPAGATOR
***********************************************************************/


void S_P11prop( int cross, double W2, double kResonance[], Matrix kresSlash, Matrix &Rprop ){
  
  double MP112 = MP11*MP11;
  double P11_Width;
  
  if( cross == 1 ){
    
    P11_Width = 0.;
    
  }else{
  
    double W = sqrt(W2);
    double EN = ( W2 - Mpi2 + MN2 )/( 2*W );
    double qcm = sqrt(EN*EN-MN2);
    
    double Width_deltapi = WdthP11*0.35;
    
    P11_Width = 3.*(fP11*fP11)* pow(W + MN,2) * (EN - MN) * qcm / (4.*Pi*Mpi2*W) + Width_deltapi;
  
  }
  

  Rprop = (1./( W2 - MP112 + I*MP11*P11_Width)) * ( kresSlash + MP11*Id );      
      
}


void P11P_current( int process, int nucleon, int decay, int Helicity, int cross, double Qsq, double W2, double Q[], double kResonance[], Matrix kpiSlash, Matrix QSlash, Matrix kresSlash, Matrix Op_R[] )
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

      double icR, icNP; 
if( cross == 0 ){
  if( process == 0 || process == 2 ){
    if( nucleon==1 ){ 
      if( decay==1 ){ icNP = sqrt(1./2.); } 
      else if( decay==2 ){ icNP = 1.; } 
    }
    else if( nucleon == 2 ){ 
      if( decay==1 ){ icNP = -sqrt(1./2.); } 
      else if( decay==2 ){ icNP = 1.; } 
     }
  }
// // // 
  else if( process == 1 ){
    if( Helicity == -1 ){
      if( nucleon==1 ){ icNP = 0.; } 
      else if( nucleon == 2 ){ 
	if( decay==1 ){ icNP = sqrt(1./2.); } 
	else if( decay==2 ){ icNP = 1.; } 
      }
   }
    else if( Helicity == 1 ){
      if( nucleon == 2 ){ icNP = 0.; } 
      else if( nucleon == 1 ){ 
	if( decay==1 ){ icNP = -sqrt(1./2.); } 
	else if( decay==2 ){ icNP = 1.; } 
      }
   }
  }
}
// // // 
else if( cross == 1 ){
  if( process == 0 || process == 2 ){
    if( nucleon==1 ){ 
      if( decay==1 ){ icNP = sqrt(1./2.); } 
      else if( decay==2 ){ icNP = 1.; } 
    }
    else if( nucleon == 2 ){ 
      if( decay==1 ){ icNP = -sqrt(1./2.); } 
      else if( decay==2 ){ icNP = 1.; } 
     }
  }
// // // 
  else if( process == 1 ){
    if( Helicity == -1 ){
      if( nucleon==1 ){ icNP = 1.; } 
      else if( nucleon == 2 ){ 
	if( decay==1 ){ icNP = -sqrt(1./2.); } 
	else if( decay==2 ){ icNP = 0.; } 
      }
    }
    else if( Helicity == 1 ){
      if( nucleon == 2 ){ icNP = 1.; } 
      else if( nucleon == 1 ){ 
	if( decay==1 ){ icNP = sqrt(1./2.); } 
	else if( decay==2 ){ icNP = 0.; } 
      }
    }
  }
}

icR = icNP;


if( icR != 0 ){
//       // -W N R- vertex      
      Matrix WNR[4];
      Gamma_WNP11( nucleon, process, decay, cross, Qsq, Q, kResonance, W2, QSlash, WNR ); 
//       // R PROPAGATOR (R)
      Matrix Rprop;
      S_P11prop( cross, W2, kResonance, kresSlash, Rprop );
//       // -R pi N- vertex 
      Matrix RNpi;
      double fact = icR * (sqrt(2.)*fP11/Mpi);
      RNpi = fact * (kpiSlash*Gamma5); 
      
      Matrix block_R;
      block_R = RNpi*Rprop;

      for( int i=0; i<4; i++ ){
	Op_R[i] = block_R*WNR[i];
      } 
  }
}


void CP11P_current( int process, int nucleon, int decay, int Helicity, int cross, double Qsq, double W2, double Q[], double kResonance[], Matrix kpiSlash, Matrix QSlash, Matrix kresSlash, Matrix Op_R_cross[] )
{  
  
// double s = pow(sMan[0],2) - pow(sMan[1],2) - pow(sMan[2],2) - pow(sMan[3],2);
  
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
      double icR, icNP; 
if( cross == 0 ){
  if( process == 0 || process == 2 ){
    if( nucleon==1 ){ 
      if( decay==1 ){ icNP = sqrt(1./2.); } 
      else if( decay==2 ){ icNP = 1.; } 
    }
    else if( nucleon == 2 ){ 
      if( decay==1 ){ icNP = -sqrt(1./2.); } 
      else if( decay==2 ){ icNP = 1.; } 
     }
  }
// // // 
  else if( process == 1 ){
    if( Helicity == -1 ){
      if( nucleon==1 ){ icNP = 0.; } 
      else if( nucleon == 2 ){ 
	if( decay==1 ){ icNP = sqrt(1./2.); } 
	else if( decay==2 ){ icNP = 1.; } 
      }
   }
    else if( Helicity == 1 ){
      if( nucleon == 2 ){ icNP = 0.; } 
      else if( nucleon == 1 ){ 
	if( decay==1 ){ icNP = -sqrt(1./2.); } 
	else if( decay==2 ){ icNP = 1.; } 
      }
   }
  }
}
// // // 
// // // 
else if( cross == 1 ){
  if( process == 0 || process == 2 ){
    if( nucleon==1 ){ 
      if( decay==1 ){ icNP = sqrt(1./2.); } 
      else if( decay==2 ){ icNP = 1.; } 
    }
    else if( nucleon == 2 ){ 
      if( decay==1 ){ icNP = -sqrt(1./2.); } 
      else if( decay==2 ){ icNP = 1.; } 
     }
  }
// // // 
  else if( process == 1 ){
    if( Helicity == -1 ){
      if( nucleon==1 ){ icNP = 1.; } 
      else if( nucleon == 2 ){ 
	if( decay==1 ){ icNP = -sqrt(1./2.); } 
	else if( decay==2 ){ icNP = 0.; } 
      }
    }
    else if( Helicity == 1 ){
      if( nucleon == 2 ){ icNP = 1.; } 
      else if( nucleon == 1 ){ 
	if( decay==1 ){ icNP = sqrt(1./2.); } 
	else if( decay==2 ){ icNP = 0.; } 
      }
    }
  }
}

icR = icNP;


if( icR != 0 ){
//       // -W N R- vertex      
      Matrix WNR_cross[4];
      Gamma_WNP11( nucleon, process, decay, cross, Qsq, Q, kResonance, W2, QSlash, WNR_cross ); 
//       // R PROPAGATOR (R)
      Matrix Rprop_cross;
      S_P11prop( cross, W2, kResonance, kresSlash, Rprop_cross );
//       // -R pi N- vertex 
      Matrix RNpi_cross;
      double fact = icR * (sqrt(2.)*fP11/Mpi); 
      RNpi_cross = fact * (kpiSlash*Gamma5); 
      
      Matrix block_R_c;
      block_R_c = Rprop_cross*RNpi_cross;
      
      for( int i=0; i<4; i++ ){
	Op_R_cross[i] = WNR_cross[i]*block_R_c;
      }
  }
  
}

#endif
