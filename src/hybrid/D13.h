#ifndef D13_H
#define D13_H

void D13_ff( double s, double u, double &D13ff)
// void D13_ff( const double &s, const double &u, double &D13ff)
{
    
  double cut_off = 800.;
  double Lam_piND = cut_off; 

  double M132=pow(MD13,2);
  double Lam_piND4 = pow( Lam_piND ,4 );
  
  double Fgauss_u = exp( -pow( u - M132, 2) / Lam_piND4 ); // Gaussian
  double Dipole_u = Lam_piND4 / ( pow( u - M132, 2) + Lam_piND4 ); // "Dipole"
  double FGaDi_u = Fgauss_u * Dipole_u;

  double Fgauss_s = exp( -pow( s - M132, 2) / Lam_piND4 ); // Gaussian
  double Dipole_s = Lam_piND4 / ( pow( s - M132, 2) + Lam_piND4 ); // "Dipole"
  double FGaDi_s = Fgauss_s * Dipole_s;

  D13ff = FGaDi_s + FGaDi_u - FGaDi_s*FGaDi_u;

}


/***************************************************************************
                          THE -W N D13- VERTEX
***************************************************************************/


// void Gamma_WND13( int nucleon, int process, int decay, int cross, const double &w, const double &q, double Q[], double kResonance[], double ki[], Matrix WND13[][4] )
void Gamma_WND13( int nucleon, int process, int decay, int cross, double Qsq, double Q[], double kResonance[], double ki[], Matrix WND13[][4] )
{

  // define the squared momentum transfer and the invariant hadronic mass
//   double Qsq = pow( q, 2 ) - pow( w, 2 );
//   double MN2 = MN*MN;


  double QsqGeV = Qsq/1.E6;
  double C3_Vp, C4_Vp, C5_Vp; 
  double C3_Vn, C4_Vn, C5_Vn;

  double DipV,M_V2;
  M_V2 = pow(0.84,2);
  DipV = pow(1.+QsqGeV/M_V2,2);
  double C5_A0;
  
// //   choose: vff=1--> Lalakulich 2006, vff=2--> Leitner 2009

  int vff = 1;
  
  if( vff == 1 ){
// // //   Lalakulich 2006 (we change the sign of proton and neutron contributions, not sure but I think it is the right way...)
    double temp=DipV*(1.+QsqGeV/(8.9*M_V2));
    
  C3_Vp = -2.95/temp;
  C4_Vp = 1.05/temp;
  C5_Vp = 0.48/DipV;

  C3_Vn = 1.13/temp;
  C4_Vn = -0.46/temp;
  C5_Vn = 0.17/DipV;
// // //   
  C5_A0 = -2.1;  //we keep the sign here as in Lalakulich
  
// //   with all this, the convention agrees with Hernandez2010 and Leitner2009
  }
  else if( vff == 2 ){ // WARNING PRODUCE SIMILAR RESULTS AT LOW Q2 (FORWARD SCATTERING ANGLES) LOT OF PROBLEMS AT HIGH Q2 !!! NOT USE THIS
// // // Leitner 2009 (the same as in Hernandez et al. PRD 87, 113009 (2013) )
    C3_Vp = -2.70/pow((1.+QsqGeV/(1.4*M_V2)),2); 
    C4_Vp = 2.62/(DipV*(1.+QsqGeV/(3.7*M_V2))); 
    C5_Vp = -1.17/(DipV*(1.+QsqGeV/(0.42*M_V2)));
  
    C3_Vn = 0.28/pow((1.+QsqGeV/(1.4*M_V2)),2); 
    C4_Vn = -1.59/(DipV*(1.+QsqGeV/(3.7*M_V2))); 
    C5_Vn = 1.96/(DipV*(1.+QsqGeV/(0.42*M_V2)));  
// // // 
    C5_A0 = -2.15;
  }
  
  double DipA,M_A2;
//   M_A2 = pow(1.0,2);
  M_A2 = 1.0;
  DipA = pow(1.+QsqGeV/M_A2,2);

  double C3_A, C4_A, C5_A, C6_A;
  C5_A = C5_A0/(DipA*(1.+QsqGeV/(3.*M_A2)));
  C3_A = 0.;
  C4_A = 0.;

  
  double C3_V, C4_V, C5_V;  
// // EM interaction  
if( process == 0 ){
  if( nucleon == 1 && decay == 1){      
      C3_V = C3_Vp; 
      C4_V = C4_Vp; 
      C5_V = C5_Vp; 
    }
  else if( nucleon == 1 && decay == 2){      
    if( cross == 0 ){
      C3_V = C3_Vp; 
      C4_V = C4_Vp; 
      C5_V = C5_Vp; 
    }
    else if( cross == 1 ){
      C3_V = C3_Vn; 
      C4_V = C4_Vn; 
      C5_V = C5_Vn; 
    }
  }    
  else if( nucleon == 2 && decay == 1){
    C3_V = C3_Vn; 
    C4_V = C4_Vn; 
    C5_V = C5_Vn; 
  }
  else if( nucleon == 2 && decay == 2){
    if( cross == 0){
      C3_V = C3_Vn; 
      C4_V = C4_Vn; 
      C5_V = C5_Vn; 
    }
    else if( cross == 1 ){
      C3_V = C3_Vp; 
      C4_V = C4_Vp; 
      C5_V = C5_Vp;
    }
  }
}
// // WNC interaction
else if( process == 2 ){
  
double C3_Vs = 0., C4_Vs = 0., C5_Vs = 0.;
double C5_As = 0., C6_As = 0.;

double wC3_Vp = 0.5*(QWeak*C3_Vp - C3_Vn - C3_Vs); 
double wC4_Vp = 0.5*(QWeak*C4_Vp - C4_Vn - C4_Vs); 
double wC5_Vp = 0.5*(QWeak*C5_Vp - C5_Vn - C5_Vs);

double wC3_Vn = 0.5*(QWeak*C3_Vn - C3_Vp - C3_Vs); 
double wC4_Vn = 0.5*(QWeak*C4_Vn - C4_Vp - C4_Vs); 
double wC5_Vn = 0.5*(QWeak*C5_Vn - C5_Vp - C5_Vs);

double wC5_Ap = 0.5*(C5_A - C5_As);
double wC5_An = 0.5*(-C5_A - C5_As);

  if( nucleon == 1 && decay == 1){      
      C3_V = wC3_Vp; 
      C4_V = wC4_Vp; 
      C5_V = wC5_Vp; 
      
      C5_A = wC5_Ap;
    }
  else if( nucleon == 1 && decay == 2){      
    if( cross == 0 ){
      C3_V = wC3_Vp; 
      C4_V = wC4_Vp; 
      C5_V = wC5_Vp; 
      
      C5_A = wC5_Ap;
    }
    else if( cross == 1 ){
      C3_V = wC3_Vn; 
      C4_V = wC4_Vn; 
      C5_V = wC5_Vn; 
      
      C5_A = wC5_An;
    }
  }    
  else if( nucleon == 2 && decay == 1){
      C3_V = wC3_Vn; 
      C4_V = wC4_Vn; 
      C5_V = wC5_Vn;
      
      C5_A = wC5_An;
  }
  else if( nucleon == 2 && decay == 2){
    if( cross == 0){
      C3_V = wC3_Vn; 
      C4_V = wC4_Vn; 
      C5_V = wC5_Vn;
      
      C5_A = wC5_An;
    }
    else if( cross == 1 ){
      C3_V = wC3_Vp; 
      C4_V = wC4_Vp; 
      C5_V = wC5_Vp;
      
      C5_A = wC5_Ap;
    }
  }
}
// // // // // //   

// // // CC interaction  
if( process == 1 ){
  C3_V = C3_Vp - C3_Vn; // C3_V(0) = -2.98
  C4_V = C4_Vp - C4_Vn; // C4_V(0) = 4.21
  C5_V = C5_Vp - C5_Vn; // C5_V(0) = -3.13 
}
// // // // // 
  
  C6_A = (pow(MN/1.E3,2))/(QsqGeV+pow(Mpi/1.E3,2))*C5_A; 

  
  
  /*
    In declaring the necessary arrays of Matrices, we distinguish between the vector part and the axial-
    vector part of the -W N D13- vertex. 
  */

  Matrix ND13_Vect[4][4];
  Matrix ND13_Ax[4][4];
  
//   Matrix q_slash;
//   q_slash = Q[0]*Gamma[0] - Q[1]*Gamma[1] - Q[2]*Gamma[2] - Q[3]*Gamma[3];
  Matrix q_slash;       //ATTENTION this quantity is not the same for delta_cross and delta, we pass -Q[] and Q[], respectively.
//   q_slash = Q[0]*Gamma[0] - Q[1]*Gamma[1] - Q[2]*Gamma[2] - Q[3]*Gamma[3];
    q_slash.M[0][0]=Q[0]       , q_slash.M[0][1]=0.         , q_slash.M[0][2]=-Q[3]        , q_slash.M[0][3]=-Q[1]+I*Q[2],
    q_slash.M[1][0]=0.         , q_slash.M[1][1]=Q[0]       , q_slash.M[1][2]=-Q[1]-I*Q[2], q_slash.M[1][3]=Q[3],
    q_slash.M[2][0]=Q[3]       , q_slash.M[2][1]=Q[1]-I*Q[2], q_slash.M[2][2]=-Q[0]       , q_slash.M[2][3]=0.,
    q_slash.M[3][0]=Q[1]+I*Q[2], q_slash.M[3][1]=-Q[3]      , q_slash.M[3][2]=0.          , q_slash.M[3][3]=-Q[0];
    
  double q_point_kD;
  q_point_kD = Q[0]*kResonance[0] - Q[1]*kResonance[1] - Q[2]*kResonance[2] - Q[3]*kResonance[3];

  double q_point_ki;
  q_point_ki = Q[0]*ki[0] - Q[1]*ki[1] - Q[2]*ki[2] - Q[3]*ki[3];
  
// // // // vector part  
  double sign;
  Matrix block3[4][4];
  double block4[4][4];
  
  for(int i=0; i<4; i++){
    if(i == 0){sign = 1.;}else{sign = -1.;}

      block3[i][i] = (C3_V/MN)  *(sign*q_slash - Q[i]*Gamma[i]);
      block4[i][i] = 1./MN2 *(sign*q_point_kD - Q[i]*kResonance[i]);
      
      ND13_Vect[i][i] = block3[i][i] + ( C4_V*block4[i][i] + C5_V*(1./MN2) * (sign*q_point_ki - Q[i]*ki[i]) )*Id;
  }

  double temp, temp1, temp2;
  for(int i=0; i<4; i++){
      temp = C3_V/MN*(-Q[i]);
      temp1 = (-Q[i])/MN2;
      temp2 = C5_V/MN2 *(-Q[i]);
    for(int j=0; j<4; j++){
       if( j != i ){
	block3[i][j] = temp *Gamma[j];
	block4[i][j] = temp1 *kResonance[j];
	
	ND13_Vect[i][j] = block3[i][j] + (C4_V*block4[i][j] + temp2*ki[j])*Id;
       }
    }
  }
// // // // // // // // // // // // //   
  
// // // // axial part 
  if( process==1 || process==2 ){
   // Assuming C3_A = C4_A = 0...
   
    double gmunu;
    double temp3;
    for(int i=0; i<4; i++){
        temp3 = C6_A/MN2*Q[i];
        for(int j=0; j<4; j++){
            
            if(i == j){gmunu = 1.;}else{gmunu = 0;}
            if(i == j && i != 0){gmunu = -1.;}
            
            ND13_Ax[i][j] = ( /*C3_A * block3[i][i] +*/  /*C4_A * block4[i][j] +*/ C5_A * gmunu + temp3*Q[j] )*Gamma5;
        }
    }  
  
  }

// // // // // // // // // 


  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){

      if(process == 0){
      WND13[i][j] = ND13_Vect[i][j] ;
      }else{
      WND13[i][j] = ( ND13_Vect[i][j] + ND13_Ax[i][j] );
      }
      
    }
  }


}


/***********************************************************************
                         THE D13 PROPAGATOR
***********************************************************************/

void S_D13( int cross, double W2, double kResonance[], Matrix kRSlash, Matrix D13[][4] )
// void S_D13( int cross, const double &W2, double kResonance[], Matrix D13[][4] )
{
//   double kResonance[4];
//   if(cross==0){
//     for( int i=0; i<4; i++ ){kResonance[i] = sMan[i];} 
//   }else{
//     for( int i=0; i<4; i++ ){kResonance[i] = uMan[i];} 
//   }
  
//   double s = pow(sMan[0],2) - pow(sMan[1],2) - pow(sMan[2],2) - pow(sMan[3],2);
//   double u = pow(uMan[0],2) - pow(uMan[1],2) - pow(uMan[2],2) - pow(uMan[3],2);
  
//   Matrix kRSlash;
//   kRSlash = kResonance[0]*Gamma[0] - kResonance[1]*Gamma[1] - kResonance[2]*Gamma[2] - kResonance[3]*Gamma[3] ;
  
  Matrix Prop_FF;
  double D13_Width;
    
//   double Mpi2 = pow(Mpi,2);
//   double MN2 = pow(MN,2);
  double MD13_2 = pow(MD13,2);
  
    if( cross == 1 ){ 
      
      D13_Width = 0.;
      
    }else{ // we compute the D13-decay width

  double W = sqrt(W2);
  double EN = ( W2 - Mpi2 + MN2 )/( 2*W );    
  double q_cm = sqrt(EN*EN-MN2); 

    double Width_deltapi = WdthD13*0.41;
    D13_Width =  ( fD13*fD13 )*( pow(q_cm,3)*( EN - MN ) )/( 4.*Pi*Mpi2*W )  + Width_deltapi;        

 }


    Prop_FF = (-1./(3.*( W2 - MD13_2 + I*MD13*D13_Width )))*( kRSlash + MD13*Id ); //the factor 1/3 is from the propagator
 
    Matrix Rarita[4][4];
//     Rarita[0][0] = ( 2. - 2./MD13_2*pow(kResonance[0],2) ) *Id ; 
//     Rarita[1][1] = ( -2. - 2./MD13_2*pow(kResonance[1],2) ) *Id ; 
//     Rarita[2][2] = ( -2. - 2./MD13_2*pow(kResonance[2],2) ) *Id ; 
//     Rarita[3][3] = ( -2. - 2./MD13_2*pow(kResonance[3],2) ) *Id ; 
    
  double a0=2. - 2./MD13_2 *pow(kResonance[0],2), a1=-2 - 2./MD13_2 *pow(kResonance[1],2), a2 = -2. - 2./MD13_2 *pow(kResonance[2],2), a3 = -2. - 2./MD13_2 *pow(kResonance[3],2);

  for(int i=0; i<4; i++){
        
    Rarita[0][0].M[i][i] = a0 ; 
    Rarita[1][1].M[i][i] = a1 ; 
    Rarita[2][2].M[i][i] = a2 ; 
    Rarita[3][3].M[i][i] = a3 ; 
      
    for(int j=0; j<4; j++){
      
      if( j != i ){
	Rarita[i][j] = mGamma_munu[i][j] - (2./MD13_2*kResonance[i]*kResonance[j])*Id + ( (kResonance[i]/MD13)*Gamma[j] - (kResonance[j]/MD13)*Gamma[i]);
      }

    }
    
  }
  
  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){

      D13[i][j] = Prop_FF*( Rarita[i][j] );
      
    }
  }

}



// void D13P_current( int process, int nucleon, int decay, int Helicity, int cross, const double &w, const double &q, const double &W2, double Q[], double sMan[], double uMan[], double ki[], double kpi[], Matrix Op_D13[] )
void D13P_current( int process, int nucleon, int decay, int Helicity, int cross, double Qsq, double W2, double Q[], double kResonance[], double ki[], double kpi[], Matrix kRSlash, Matrix Op_D13[] )
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
      double icNP, icD13;
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

icD13 = icNP;



if( icD13 == 0 ){
   for( int i=0; i<4; i++ ){
    Op_D13[i] = 0*Id;
   }  
}
else{
//       // -W N D13- vertex 
      Matrix WND13[4][4]; 
      Gamma_WND13( nucleon, process, decay, cross, Qsq, Q, kResonance, ki, WND13 );
//       // D13 propagator
      Matrix D13[4][4];
      S_D13( cross, W2, kResonance, kRSlash, D13 ); 
// //  D13[i][j] has been defined in CONTRAVARIANT (upper indices) notation but we need it in covariant (down indices)
      for( int i=1; i<4; i++ ){
	  D13[0][i] = (-1.)*D13[0][i];
	  D13[i][0] = (-1.)*D13[i][0];
	}      
//       // -D13 pi N- vertex 
      Matrix D13Npi[4];
      Matrix fact;
      fact = (icD13*(fD13*sqrt(2.)/Mpi))*Gamma5;
      for( int i=0; i<4; i++ ){
	D13Npi[i] =  kpi[i]*fact;
      }
      
// // // 

//       for( int i=0; i<4; i++ ){
// 	for( int k=0; k<4; k++ ){
// 	  for( int l=0; l<4; l++ ){
// 	    Op_D13[i] = Op_D13[i] + D13Npi[k]*D13[k][l]*WND13[l][i];
// 	  }
// 	}
//       }

    Matrix temp;
    for( int k=0; k<4; k++ ){
        for( int l=0; l<4; l++ ){
            temp = D13Npi[k]*D13[k][l];
            for( int i=0; i<4; i++ ){
          
                Op_D13[i] = Op_D13[i] + temp*WND13[l][i];
            
            }
        }
    }

      
  }
}



// void CD13P_current( int process, int nucleon, int decay, int Helicity, int cross, const double &w, const double &q, const double &W_cross2, double minusQ[], double sMan[], double uMan[], double kN[], double kpi[], Matrix Op_D13_cross[] )
void CD13P_current( int process, int nucleon, int decay, int Helicity, int cross, double Qsq, double W2, double Q[], double kResonance[], double kN[], double kpi[], Matrix kRSlash, Matrix Op_D13_cross[] )
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
      double icNP, icD13;
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

icD13 = icNP;


if( icD13 == 0 ){
  for( int i=0; i<4; i++ ){
    Op_D13_cross[i] = 0.*Id;
  }
}else{
//       // -W N D13- vertex 
      Matrix WND13_cross[4][4]; 
      Gamma_WND13( nucleon, process, decay, cross, Qsq, Q, kResonance, kN, WND13_cross );
      
//       for( int i=0; i<4; i++ ){ for( int j=0; j<4; j++ ){ cout << WND13_cross[i][j] << endl; } }
// // D13 propagator
      Matrix D13_cross[4][4];
      S_D13( cross, W2, kResonance, kRSlash, D13_cross );


      for( int i=1; i<4; i++ ){
	  D13_cross[0][i] = (-1.)*D13_cross[0][i];
	  D13_cross[i][0] = (-1.)*D13_cross[i][0];
	}     
//       // -D13 pi N- vertex 
      Matrix D13Npi_cross[4];
      Matrix fact;
      fact =  (icD13*(-fD13*sqrt(2.)/Mpi))*Gamma5;
      for( int i=0; i<4; i++ ){
	D13Npi_cross[i] = kpi[i] * fact;
      }
 
//       Matrix block_D13_c[4][4];
//       for( int i=0; i<4; i++ ){
// 	for( int j=0; j<4; j++ ){
// 	  block_D13_c[i][j] = Gamma[0] * ConjMatrix( WND13_cross[j][i] ) * Gamma[0];
// 	}
//       }
      Matrix block_D13_c[4][4];
      for( int i=0; i<4; i++ ){
	for( int j=0; j<4; j++ ){
	  block_D13_c[i][j] = ConjMatrix( WND13_cross[j][i] );
          
                block_D13_c[i][j].M[0][2] = -block_D13_c[i][j].M[0][2];
                block_D13_c[i][j].M[0][3] = -block_D13_c[i][j].M[0][3];
                block_D13_c[i][j].M[1][2] = -block_D13_c[i][j].M[1][2];
                block_D13_c[i][j].M[1][3] = -block_D13_c[i][j].M[1][3];
                
                block_D13_c[i][j].M[2][0] = -block_D13_c[i][j].M[2][0];
                block_D13_c[i][j].M[2][1] = -block_D13_c[i][j].M[2][1];
                block_D13_c[i][j].M[3][0] = -block_D13_c[i][j].M[3][0];
                block_D13_c[i][j].M[3][1] = -block_D13_c[i][j].M[3][1];
                
	}
      }



//       for( int i=0; i<4; i++ ){
// 	for( int k=0; k<4; k++ ){
// 	  for( int l=0; l<4; l++ ){
// 	    Op_D13_cross[i] = Op_D13_cross[i] + block_D13_c[i][l]*D13_cross[l][k]*D13Npi_cross[k];
// 	  }
// 	}
//       }
      
    Matrix temp;
    for( int k=0; k<4; k++ ){
        for( int l=0; l<4; l++ ){
            temp = D13_cross[l][k]*D13Npi_cross[k];
            for( int i=0; i<4; i++ ){
                
                Op_D13_cross[i] = Op_D13_cross[i] + block_D13_c[i][l]*temp;
                
            }
        }
    }
      
      
  }
  
}

#endif