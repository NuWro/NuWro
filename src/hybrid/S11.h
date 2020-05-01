#ifndef S11_H
#define S11_H

// /***************************************************************************
//                           THE -W N R- VERTEX
// ***************************************************************************/

void S11_ff(double s, double u, double &S11ff)
{

  double cut_off = 1200.;
  double Lam_piND = cut_off; 
  
  double Lam_piND4=pow( Lam_piND ,4 );
  double MS112 = pow(MS11,2);

  double Fgauss_u = exp( -pow( u - MS112, 2) / Lam_piND4 ); // Gaussian
  double Dipole_u = Lam_piND4 / ( pow( u - MS112, 2) + Lam_piND4 ); // "Dipole"
  double FGaDi_u = Fgauss_u * Dipole_u;

  double Fgauss_s = exp( -pow( s - MS112, 2) / Lam_piND4 ); // Gaussian
  double Dipole_s = Lam_piND4 / ( pow( s - MS112, 2) + Lam_piND4 ); // "Dipole"
  double FGaDi_s = Fgauss_s * Dipole_s;

  S11ff = FGaDi_s + FGaDi_u - FGaDi_s*FGaDi_u;
  
}


void Gamma_WNS11( int nucleon, int process, int decay, int cross, double Qsq, double Q[], double kResonance[], Matrix QSlash, Matrix WNR[] )
{ 
  
  double QsqGeV, xmu, tau;
  QsqGeV = Qsq/1.E6;
  xmu = MN + MS11;
  
//   double W = sqrt( pow(kResonance[0],2) - pow(kResonance[1],2) - pow(kResonance[2],2) - pow(kResonance[3],2) );
  
  double DipV, DipA, MA2, MV2;
  MA2 = 1.1025; //MA=1.05 GeV
  MV2 = 0.7100; //MV=0.84 GeV
  DipV = pow( 1+QsqGeV/MV2, 2 ); 
  DipA = pow( 1+QsqGeV/MA2, 2 );  
    
  
  double F1p, F2p, F1n, F2n, GA0, GP0;  
  
  int Vff = 1;
  
  if( Vff == 1 ){
// // // Lalakulich 2006 WARNING we consider a relative sign in the proton form factor in order to use our convention of isovector form factors   
  F1p = -2.0/(DipV*(1.+ QsqGeV/(1.2*MV2))) * ( 1.+7.2*log(1.+QsqGeV/1.0) );
  F2p = -0.84/DipV*( 1.+ 0.11*log(1.+QsqGeV/1.0) );
  F1n = -F1p;
  F2n = -F2p;

  GA0 = -0.21 / ( DipA * (1. + QsqGeV/(3.*MA2)) ) ;
// // // // // // // // // // // //   
  }
  else if( Vff == 2 ){ 
// // // Leitner 2009, WARNING! In Leitner et al. 2009 they did not published the Qsq-dependence of the vector form factors, so this is only a test: dont use it for calculations!!!
  F1p = 0.85/DipV;
  F2p = 0.46/DipV;
  F1n = -0.02/DipV;
  F2n = -0.35/DipV;
  
    GA0 = -0.23 / DipA;
// // // // // // // // // // // //   
  }
  
  
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
else if( process == 2 ){
    
  double F1s = 0., F2s = 0.;
  double GAs =0.;
  
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
else if( process == 1 ){
  F1 = (F1p - F1n);
  F2 = (F2p - F2n);
}
// // // // // //   

  GP0 = (MS11-MN)/(Qsq + Mpi2) * GA0; //MeV-1

  
  
//   Matrix QSlash; 
// //   QSlash = Q[0]*Gamma[0] - Q[1]*Gamma[1] - Q[2]*Gamma[2] - Q[3]*Gamma[3];
// 
// QSlash.M[0][0]=Q[0]       , QSlash.M[0][1]=0.         , QSlash.M[0][2]=-Q[3]        , QSlash.M[0][3]=-Q[1]+I*Q[2],
// QSlash.M[1][0]=0.         , QSlash.M[1][1]=Q[0]       , QSlash.M[1][2]=-Q[1]-I*Q[2], QSlash.M[1][3]=Q[3],
// QSlash.M[2][0]=Q[3]       , QSlash.M[2][1]=Q[1]-I*Q[2], QSlash.M[2][2]=-Q[0]       , QSlash.M[2][3]=0.,
// QSlash.M[3][0]=Q[1]+I*Q[2], QSlash.M[3][1]=-Q[3]      , QSlash.M[3][2]=0.          , QSlash.M[3][3]=-Q[0];  
  
//   Matrix WNR_V[4];
//   Matrix WNR_A[4];  
// // // // // // // Vector part // // // // // // // // // // // // // // 
//   for(int i=0; i<4; i++){
//     WNR_V[i] = F1/(pow(xmu,2)) * (Qsq*Gamma[i] + Q[i]*QSlash) - F2/(xmu)*(Gamma[i]*QSlash-QSlash*Gamma[i]);
//   }
//   
//   
// // // // // // // Axial part // // // // // // // // // // // // // // 
// 
//   if( process == 1 || process == 2 ){
//     for(int i=0; i<4; i++){
//       WNR_A[i] = GA0*Gamma_mu5[i] + (GP0* Q[i])*Gamma5;
//     }
//   }
// // // // // // // // // // // // // // // // // // // // // // // // // // 
//   for(int i=0; i<4; i++){
// 
//     if( process == 0 ){
//       WNR[i] = WNR_V[i]*Gamma5 ;
//     }
//     if( process == 1 || process == 2 ){    
//       WNR[i] = ( WNR_V[i] - WNR_A[i] )*Gamma5;
//     }
//     
//   }  
     



    if( process == 0 ){
        
  for(int i=0; i<4; i++){
      WNR[i] = ( F1/(pow(xmu,2)) * (Qsq*Gamma[i] + Q[i]*QSlash) - F2/(xmu)*(Gamma[i]*QSlash-QSlash*Gamma[i]) )*Gamma5 ;
  }
  
    }else{
        
  for(int i=0; i<4; i++){
      WNR[i] = ( F1/(pow(xmu,2)) * (Qsq*Gamma[i] + Q[i]*QSlash) - F2/(xmu)*(Gamma[i]*QSlash-QSlash*Gamma[i]) - ( GA0*Gamma_mu5[i] + (GP0* Q[i])*Gamma5 ) )*Gamma5;
  }
  
    }
    

}


/***********************************************************************
                         THE S11 PROPAGATOR
***********************************************************************/


void S_S11prop( int cross, double W2, double kResonance[], Matrix kresSlash, Matrix &Rprop ){
  
//   double kResonance[4];
//   if(cross==0){
//     for( int i=0; i<4; i++ ){kResonance[i] = sMan[i];} 
//   }else{
//     for( int i=0; i<4; i++ ){kResonance[i] = uMan[i];} 
//   }
//   
//   double s = pow(sMan[0],2) - pow(sMan[1],2) - pow(sMan[2],2) - pow(sMan[3],2);
//   double u = pow(uMan[0],2) - pow(uMan[1],2) - pow(uMan[2],2) - pow(uMan[3],2);
//   double W = sqrt(s);
  
//   double Mpi2 = Mpi*Mpi;
//   double MN2 = MN*MN;
  double MS112 = MS11*MS11;
  double S11_Width;
  
  if( cross == 1 ){
    
    S11_Width = 0.;
    
  }else{
  
    double W = sqrt(W2);
    double EN = ( W2 - Mpi2 + MN2 )/( 2*W );
    double qcm = sqrt(EN*EN-MN2);
    
    double Width_deltapi = WdthS11*0.55;
    
    S11_Width = 3.*(fS11*fS11)* pow(W - MN,2) * (EN + MN) * qcm / (4.*Pi*Mpi2*W) + Width_deltapi;
  
  }

//   Rprop = (1./( W2 - MS112 + I*MS11*S11_Width)) *( kResonance[0]*Gamma[0] - kResonance[1]*Gamma[1] - kResonance[2]*Gamma[2] - kResonance[3]*Gamma[3] + MS11*Id );
  Rprop = (1./( W2 - MS112 + I*MS11*S11_Width)) * ( kresSlash + MS11*Id );         
      
}


void S11P_current( int process, int nucleon, int decay, int Helicity, int cross, double Qsq, double W2, double Q[], double kResonance[], Matrix kpiSlash, Matrix QSlash, Matrix kresSlash, Matrix Op_R[] )
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


// if( icR == 0 ){
//   for( int i=0; i<4; i++ ){
//     Op_R[i] = 0*Id;
//   }
// }
// else{
if( icR != 0 ){
//       // -W N R- vertex      
      Matrix WNR[4];
      Gamma_WNS11( nucleon, process, decay, cross, Qsq, Q, kResonance, QSlash, WNR ); 
//       // R PROPAGATOR (R)
      Matrix Rprop;
      S_S11prop( cross, W2, kResonance, kresSlash, Rprop );
//       // -R pi N- vertex 
      Matrix RNpi;
      double fact = icR * (sqrt(2.)*fS11/Mpi);
      RNpi =  fact * kpiSlash; 
      
      Matrix block_R;
      block_R = RNpi*Rprop;

      for( int i=0; i<4; i++ ){
	Op_R[i] = block_R*WNR[i];
      } 
  }
}


void CS11P_current( int process, int nucleon, int decay, int Helicity, int cross, double Qsq, double W2, double Q[], double kResonance[], Matrix kpiSlash, Matrix QSlash, Matrix kresSlash, Matrix Op_R_cross[] )
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
      if( nucleon==1 ){ icNP = 0.; } 
      if( nucleon == 2 ){ 
	if( decay==1 ){ icNP = sqrt(1./2.); } 
	if( decay==2 ){ icNP = 1.; } 
      }
   }
    if( Helicity == 1 ){
      if( nucleon == 2 ){ icNP = 0.; } 
      if( nucleon == 1 ){ 
	if( decay==1 ){ icNP = -sqrt(1./2.); } 
	if( decay==2 ){ icNP = 1.; } 
      }
   }
  }
}
// // // 
// // // 
if( cross == 1 ){
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
      if( nucleon==1 ){ icNP = 1.; } 
      if( nucleon == 2 ){ 
	if( decay==1 ){ icNP = -sqrt(1./2.); } 
	if( decay==2 ){ icNP = 0.; } 
      }
    }
    if( Helicity == 1 ){
      if( nucleon == 2 ){ icNP = 1.; } 
      if( nucleon == 1 ){ 
	if( decay==1 ){ icNP = sqrt(1./2.); } 
	if( decay==2 ){ icNP = 0.; } 
      }
    }
  }
}

icR = icNP;


// if( icR == 0 ){
//   for( int i=0; i<4; i++ ){
//     Op_R_cross[i] = 0*Id;
//   } 
// }
// else{
if( icR != 0){
//       // -W N R- vertex      
      Matrix WNR_cross[4];
      Gamma_WNS11( nucleon, process, decay, cross, Qsq, Q, kResonance, QSlash, WNR_cross ); 
//       // R PROPAGATOR (R)
      Matrix Rprop_cross;
      S_S11prop( cross, W2, kResonance, kresSlash, Rprop_cross );
//       // -R pi N- vertex 
      Matrix RNpi_cross;     
      double fact = icR * (sqrt(2.)*fS11/Mpi) ;
      RNpi_cross = fact * kpiSlash; 
      
      Matrix block_R_c;
      block_R_c = Rprop_cross*RNpi_cross;
      
      for( int i=0; i<4; i++ ){
	Op_R_cross[i] = WNR_cross[i]*block_R_c;
      }
  }
  
}

#endif