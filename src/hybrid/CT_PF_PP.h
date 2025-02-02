#ifndef CT_PF_PP_H
#define CT_PF_PP_H


void formfactor(double Qsq, double t, double &F1V, double &F2V, double &F1p, double &F1n, double &F2p, double &F2n, double &F_rho)
{
    
  double QsqGeV, xmu, tau;
  QsqGeV = Qsq/1.E6;
  xmu = 2.*MN;
  tau = Qsq/(xmu*xmu);
    
  double DipV, DipA, MA2, MV2;
  MA2 = 1.1025; //MA=1.05 GeV
  MV2 = 0.7100; //MV=0.84 GeV
  DipV = pow( 1.+QsqGeV/MV2, 2 ); 
  DipA = pow( 1.+QsqGeV/MA2, 2 );  
    
  
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
    
    F1p = (GEp+tau*GMp)/(1.+tau);
    F2p = (GMp-GEp)/(1.+tau);
    F1n = (GEn+tau*GMn)/(1.+tau);
    F2n = (GMn-GEn)/(1.+tau);

    F1V = (F1p - F1n);
    F2V = (F2p - F2n);

// // // // // // // //     
    
    double mrho = 775.8; // MeV    
    
    F_rho = 1./( 1. - t/(mrho*mrho) );

}


void CT_PF_PP( int process, int nucleon, int decay, int Helicity, double Q[], double kpi[], double t, double Qsq, Matrix QSlash, Matrix Op_CT[], Matrix Op_PF[], Matrix Op_PP[] )
{
  
//double Mpi2 = Mpi*Mpi, MN2 = MN*MN;

double F1V, F2V, F_rho, F1p, F1n, F2p, F2n;
      formfactor(Qsq, t, F1V, F2V, F1p, F1n, F2p, F2n, F_rho);
      
double F_CT, F_PF;

if( process == 2 ){
  F1V = (1.-2.*Sin2W)*F1V;
}

F_CT = F1V;
F_PF = F_CT;
// // // // // // // // // // // //   


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
    
double icCT, icPP, icPF;
if( process == 1 ){
  if(Helicity == -1){
    if( nucleon == 1 ){ icCT = 1.;}
    if( nucleon == 2 ){ 
      if( decay == 1 ){ icCT = -sqrt(2.); }
      if( decay == 2 ){ icCT = -1.; } 
    }
  }
  if(Helicity == 1){
    if( nucleon == 2 ){ icCT = 1.;}
    if( nucleon == 1 ){ 
      if( decay == 1 ){ icCT = sqrt(2.); }
      if( decay == 2 ){ icCT = -1.; } 
    }
  }
}
if( process == 0 || process == 2 ){
  if( nucleon == 1 ){
    if( decay == 1 ){ icCT = 0;}
    if( decay == 2 ){ icCT = -1.;}
  }
  if( nucleon == 2 ){
    if( decay == 1 ){ icCT = 0;}
    if( decay == 2 ){ icCT = 1.;}
  }  
}
    
  icPP = icCT;
  icPF = icCT;


// // //   CONTACT TERM 
double fact_ct = - icCT/(sqrt(2.)*fpi);

if( icCT == 0 ){
  for( int i=0; i<4; i++ ){
    Op_CT[i] = 0.*Id;
  }
}else{
    for( int i=0; i<4; i++ ){
      for (int l=0; l < 4;l++)
      {
        for (int k=0;k<4;k++)
        {
          Op_CT[i].M[l][k] = (fact_ct*gA*F_CT) * (Gamma_mu5[i].M[l][k]);
          if( process == 1 || process == 2 ){
              Op_CT[i].M[l][k] -= (fact_ct*F_rho) * (Gamma[i].M[l][k]);
          }
        }
      }
    }
}
// // // // // // // // // // // // // // 

if( process == 0 ){icPP=0;} 
// 
if( icPP == 0 ){
  for( int i=0; i<4; i++ ){
    Op_PP[i] = 0.*Id;
  }
}else{
  double fact_pp;
  fact_pp = (-icPP/(sqrt(2.)*fpi)*F_rho/(-Qsq-Mpi2));
  for( int i=0; i<4; i++ ){
    for( int l=0 ; l<4 ; l++)
    {
        for( int k=0 ; k < 4; k++)
        {
            Op_PP[i].M[l][k] = Q[i] * fact_pp * (QSlash.M[l][k]);
        }
    }
  }
}
// // // // // // // // // // // // // // 
// 
if( icPF == 0 ){
      for( int i=0; i<4; i++ ){
	Op_PF[i] = 0.*Id;
      }
}else{
  double fact_pf = -icPF*gA/(sqrt(2.)*fpi)*F_PF / (t-Mpi2);
  
    complex<double> tslashmu5 = 0;
    for (int l=0;l<4;l++)
    {
        for(int k=0; k<4; k++)
        {
            tslashmu5= (Q[0]-kpi[0])*Gamma_mu5[0].M[l][k] - (Q[1]-kpi[1])*Gamma_mu5[1].M[l][k] - (Q[2]-kpi[2])*Gamma_mu5[2].M[l][k] - (Q[3]-kpi[3])*Gamma_mu5[3].M[l][k];
            for( int i=0; i<4; i++ ){
            Op_PF[i].M[l][k] = (fact_pf *(2.*kpi[i]-Q[i])) *tslashmu5;

            }
        }
    }
  }

}

#endif
