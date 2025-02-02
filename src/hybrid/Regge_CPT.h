#ifndef REGGE_H
#define REGGE_H


void Pion_Nucleon_FF(int process, int nucleon, int decay, int ChPT, double Qsq, double t, double s, double u, complex<double> &Fbpipi, complex<double> &wFCT, complex<double> &F1, complex<double> &F1_c, complex<double> &GA0, complex<double> &GA0_c, complex<double> &PRegge_NP, complex<double> &PRegge_rho)
{
  double QsqGeV = Qsq/1.E6;
//   double Mpi2 = Mpi*Mpi;
//   double MN2 = MN*MN;
  
// // // // // Pion in flight (PF) 
  
  double Lamb_bpipi=655.;

//   double Lamb_bpipi;  
//   if( QsqGeV < 0.4 ){ Lamb_bpipi = 775.; }
//   if( QsqGeV < 1.5 && QsqGeV > 0.4 ){ Lamb_bpipi = 630.; }
//   if( QsqGeV > 1.5 ){ Lamb_bpipi = 680.; }
  
  double fbpipi = 1./(1.+Qsq/pow(Lamb_bpipi,2));
  double alpha_prime = 0.74E-6;// MeV^-2 ; 0.74 GeV^-2 
  double alpha_pi = alpha_prime*(t-Mpi2);
  double m_alpha_pi = -alpha_pi;
  
// // // PHASE // // //
  complex<double> Rphase_PF;  
  double param;   
    
// // // VR choice of the phase for electroproduction
//   if( (process == 0 || process == 2) && nucleon == 1 && decay == 2 ){ //choice for electroproduction of p(e,e'pi+)n
//     param = 0., Rphase_PF =  param + (1.-param)*(cos(Pi*alpha_pi) - I*sin(Pi*alpha_pi)) ;     
//   } 
//   if( (process == 0 || process == 2) && nucleon == 2 && decay == 2 ){ //choice for electroproduction of n(e,e'pi-)p
//     param = 1., Rphase_PF =  param + (1.-param)*(cos(Pi*alpha_pi) - I*sin(Pi*alpha_pi)) ; 
//   } 
//   if( (process == 0 || process == 2) && decay == 1 ){ //choice for electroproduction of n(p)(e,e'pi0)n(p)
//     param = 0., Rphase_PF =  param + (1.-param)*(cos(Pi*alpha_pi) - I*sin(Pi*alpha_pi)) ; // Pion-exchanged term does not contribute 
//   } 
// // // // 

  
  
// // // Resonances (NP and CNP)

// // // // Vrancx&Ryckebusch form factors // // // // 
  double MV = 840.; //MeV
  double MA = 1050.; //MeV
  double Minf = 2194.; //MeV
  double Minf_A = 7200.; //MeV
  double LambdaV, LambdaA, LambdaV_c, LambdaA_c;
  
    LambdaV = MV + (Minf-MV)*(1. - MN2/s); 
    LambdaA = MA + (Minf_A-MA)*(1. - MN2/s);
    
    LambdaV_c = MV + (Minf-MV)*(1. - MN2/(2.*MN2-u) ); 
    LambdaA_c = MA + (Minf_A-MA)*(1. - MN2/(2.*MN2-u) );
   
  complex<double> F1p_VR = 1./pow( 1. + Qsq/pow(LambdaV,2) ,2);
  complex<double> GAp_VR = gA/pow( 1. + Qsq/pow(LambdaA,2) ,2);
// 
  complex<double> F1p_c_VR = 1./pow( 1. + Qsq/pow(LambdaV_c,2) ,2);
  complex<double> GAp_c_VR = gA/pow( 1. + Qsq/pow(LambdaA_c,2) ,2);

// 
// // // // // Kaskulov&Mosel form factors // // // // 
//   
//   complex<double> F1p_KM, F1p_c_KM, GAp_KM, GAp_c_KM;
//   double xi = 0.4;
//   
//   F1p_KM = ( s * log( xi*Qsq/MN2 + 1. ) * (2.*xi*Qsq+s)/pow(xi*Qsq,2) - s*(xi*Qsq+s)/(xi*Qsq*(xi*Qsq+MN2)) + log((s-MN2)/MN2) - I*Pi ) / ( pow(xi*Qsq/s + 1.,2) * ( (s*s + 2.*s*MN2)/(2.*MN2*MN2) + log( (s-MN2)/MN2 ) - I*Pi ) );
//   
//   GAp_KM = F1p_KM * gA;
// 
//   F1p_c_KM = ( u * log( xi*Qsq/MN2 + 1. ) * (2.*xi*Qsq+u)/pow(xi*Qsq,2) - u*(xi*Qsq+u)/(xi*Qsq*(xi*Qsq+MN2)) + log((MN2-u)/MN2) ) / ( pow(xi*Qsq/u + 1.,2) * ( (u*u + 2.*u*MN2)/(2.*MN2*MN2) + log( (MN2-u)/MN2 ) ) );
// 
//   GAp_c_KM = F1p_c_KM * gA;
// // // // // // // // // // // // // // // // // //

  complex<double> F1p, F1p_c, GAp, GAp_c;
//   F1p = F1p_KM, F1p_c = F1p_c_KM, GAp = GAp_KM, GAp_c = GAp_c_KM; //Kask&Mosel FF
  F1p = F1p_VR, F1p_c = F1p_c_VR, GAp = GAp_VR, GAp_c = GAp_c_VR; //Vrancx&Ryckebusch FF  
// // 
  
  double aa;
  if(ChPT == 1){aa=0.;}else{aa =2.4;}
  
  double alpha_prime_NP = alpha_prime/(1.+aa*Qsq/s);// MeV^-2 ; 0.74 GeV^-2 
  double alpha_pi_NP = alpha_prime_NP*(t-Mpi2);
  double m_alpha_pi_NP = -alpha_pi_NP;  

  complex<double> Rphase_NP, Rphase_NP_c;
//   if( (process == 0 || process == 2) && decay != 1){
//     param = 0., Rphase_NP =  param + (1.-param)*(cos(Pi*alpha_pi_NP) - I*sin(Pi*alpha_pi_NP)) ;
//     param = 1., Rphase_NP_c =  param + (1.-param)*(cos(Pi*alpha_pi_NP) - I*sin(Pi*alpha_pi_NP)) ;
//   }
// 
//   if( (process == 0 || process == 2) && decay == 1 ){
//     param = 1., Rphase_NP =  param + (1.-param)*(cos(Pi*alpha_pi_NP) - I*sin(Pi*alpha_pi_NP)) ;
//     param = 1., Rphase_NP_c =  param + (1.-param)*(cos(Pi*alpha_pi_NP) - I*sin(Pi*alpha_pi_NP)) ;
//   }
//   
//   if( process == 1 ){ //CHECK THIS, I did A COMPLETELY ARBITRARY CHOICE
//     param = 1., Rphase_PF =  param + (1.-param)*(cos(Pi*alpha_pi) - I*sin(Pi*alpha_pi)) ;
//     param = 1., Rphase_NP =  param + (1.-param)*(cos(Pi*alpha_pi_NP) - I*sin(Pi*alpha_pi_NP)) ;
//     Rphase_NP_c = Rphase_NP;
//   }
  
 
  //ALEXIS: set rotating phase !
// // // Constant phase for all terms
    param = 1., Rphase_PF =  param + (1.-param)*(cos(Pi*alpha_pi) - I*sin(Pi*alpha_pi)) ;
    param = 1., Rphase_NP =  param + (1.-param)*(cos(Pi*alpha_pi_NP) - I*sin(Pi*alpha_pi_NP)) ;
    Rphase_NP_c = Rphase_NP;
// // // // // // // // //  
 
// // // // // Pion in flight (PF) 
  complex<double> PRegge_PF;

  PRegge_PF = -alpha_prime * Rphase_PF * tgamma(m_alpha_pi) * pow(alpha_prime*s,alpha_pi);
  
  PRegge_PF = PRegge_PF * (t-Mpi2);


  Fbpipi = fbpipi * PRegge_PF;  
 
 
// // // Resonances (NP and CNP) 
//   complex<double> PRegge_NP;
  PRegge_NP = (-alpha_prime_NP) * tgamma(m_alpha_pi_NP) * pow(alpha_prime_NP*s,alpha_pi_NP)* (t-Mpi2); 
  
  complex<double> Fp, Fp_c, GA, GA_c;
  Fp = F1p * PRegge_NP * Rphase_NP ;
//   GA = GAp * PRegge_NP * Rphase_NP ;
  
  Fp_c = F1p_c * PRegge_NP * Rphase_NP_c ;
//   GA_c = GAp_c * PRegge_NP * Rphase_NP_c ;

  
  
// // // // // //  rho propagator for axial current // // // // // //
  double mrho2 = pow(775.8,2);
  double alpha_prime_rho = 0.85E-6;
  double alpha_rho = 0.53 + alpha_prime_rho * t;
  complex<double> phase_rho = 1.;
//Alexis: set phase to rotating!
    param = 1., phase_rho =  param + (1.-param)*(cos(Pi*alpha_rho) - I*sin(Pi*alpha_rho)) ;
  PRegge_rho = -alpha_prime_rho * phase_rho * tgamma(1.-alpha_rho) * pow(alpha_prime_rho*s,alpha_rho-1.) * (t - mrho2);
  
  GA = GAp * PRegge_rho ;  
  GA_c = GAp_c * PRegge_rho ;
// // // // // // // // // // // // // // // // // // // // // // //   
  
  
  
  if( process == 0 ){
    if( nucleon == 1 && decay == 1){ //p(e,e'pi0)p   
//       F1 = Fp;
//       F1_c = Fp_c;
      F1 = 0.;
      F1_c = 0.;
      
    }else if( nucleon == 1 && decay == 2){ //p(e,e'pi+)n     
      F1 = Fp;
      F1_c = 0.;
      
    }else if( nucleon == 2 && decay == 1){ //n(e,e'pi0)n
      F1 = 0.;
      F1_c = 0.;
      
    }else if( nucleon == 2 && decay == 2){ //n(e,e'pi-)p
      F1 = 0.;
      F1_c = Fp_c;
    }
  }
  
if( process == 2 ){
  
  complex<double> wF1p = 0.5*(QWeak*Fp); 
  complex<double> wF1n = 0.5*(-Fp); 
  complex<double> wGAp = 0.5*(GA);
  complex<double> wGAn = 0.5*(-GA);

  complex<double> wF1p_c = 0.5*(QWeak*Fp_c); 
  complex<double> wF1n_c = 0.5*(-Fp_c); 
  complex<double> wGAp_c = 0.5*(GA_c);
  complex<double> wGAn_c = 0.5*(-GA_c);
  
  if( nucleon == 1 && decay == 1){ //p(e,e'pi0)p
    F1 = wF1p;
    F1_c = wF1p_c;

    GA0 = wGAp;
    GA0_c = wGAp_c;
  }
  if( nucleon == 1 && decay == 2){ //p(e,e'pi+)n
    F1 = wF1p;
    F1_c = wF1n_c;

    wFCT = wF1p-wF1n_c;

    GA0 = wGAp;
    GA0_c = wGAn_c;
  }
  if( nucleon == 2 && decay == 1){ //n(e,e'pi0)n
    F1 = wF1n;
    F1_c = wF1n_c;
    
    GA0 = wGAn;
    GA0_c = wGAn_c;
  }
  if( nucleon == 2 && decay == 2){ //n(e,e'pi-)p
    F1 = wF1n;
    F1_c = wF1p_c;

    wFCT = wF1p_c-wF1n;

    GA0 = wGAn;
    GA0_c = wGAp_c;
  }
  
  Fbpipi = (1.-2.*Sin2W)*Fbpipi;

}
// // // //     
    
  
// // 
  if( process == 1 ){
    F1 = Fp;
    F1_c = Fp_c;
    
    GA0 = GA;
    GA0_c = GA_c;
  }
// // // 
  
  
}


void Regge( int process, int nucleon, int decay, int Helicity, double Qsq, double P[], double pN[], double Q[], double kpi[], Matrix QSlash, Matrix kpiSlash, double u, double s, double t, Matrix sSlash, Matrix uSlash, Matrix tSlash, Matrix PF_regge[], Matrix NPv_regge[], Matrix NPa_regge[] , Matrix CNPv_regge[], Matrix CNPa_regge[] , Matrix CTa_regge[], Matrix CTv_regge[], Matrix PP_regge[] )
{  
  
  double Mpi2 = Mpi*Mpi;
  double MN2 = MN*MN;   

// // 

int ChPT = 1; //if ChPT = 1 then we use Chiral Perturbation Theory operators (PF + CT + NP(CNP))
// // 
    
  complex<double> Fbpipi, wFCT, Fp, Fp_c, GA, GA_c;
  complex<double> PRegge_NP, PRegge_rho;
    Pion_Nucleon_FF(process, nucleon, decay, ChPT, Qsq, t, s, u, Fbpipi, wFCT, Fp, Fp_c, GA, GA_c, PRegge_NP, PRegge_rho );
    
// // // // // // 

    
// double gPiNN = 13.4; //Kask&Mosel choice
double gPiNN = 13.; // VR choice  // not used if using Chiral model

// // Hadronic form factor VR model (better agreement at high t (it is better to do not include this for -t<0.5 ) )
double Fmpi=1.;
//   if( (-t)/1.E6 > 0.5 ){
//     double Lambda_mpi=800.;
//     Fmpi = 1./(1.+ (Mpi2-t)/pow(Lambda_mpi,2));
//     gPiNN = gPiNN*Fmpi;
//   }
  
  Fbpipi = Fbpipi*Fmpi, wFCT = wFCT*Fmpi, Fp=Fp*Fmpi, Fp_c=Fp_c*Fmpi, GA=GA*Fmpi, GA_c = GA_c*Fmpi;
    
// // // // // // 
  
  
// // // // // //   

double icPF;
if( process == 1 ){
  if(Helicity == -1){
    if( nucleon == 1 ){ icPF = 1.;}
    if( nucleon == 2 ){ 
      if( decay == 1 ){ icPF = -sqrt(2.); }
      if( decay == 2 ){ icPF = -1.; } 
    }
  }
    
//   HELICITY == 1 --> W^- induced 1-pion production 
// 	nucleon = 2 --> neutron initial state
// 		decay = 1 or 2 --> n + pi^-
// 	nucleon = 1 --> proton initial state
// 		decay = 1 --> n + pi^0
// 		decay = 2 --> p + pi^+

  if(Helicity == 1){
    if( nucleon == 2 ){ icPF = 1.;} //nub + n -->n + pi^-
    if( nucleon == 1 ){ 
      if( decay == 1 ){ icPF = sqrt(2.); } //nub + p --> n + pi^0
      if( decay == 2 ){ icPF = -1.; } //nub + p --> p + pi^-
    }
  }
}

if( process == 0 || process == 2 ){
  if( nucleon == 1 ){
    if( decay == 1 ){ icPF = 0;}
    if( decay == 2 ){ icPF = -1.;}
  }
  if( nucleon == 2 ){
    if( decay == 1 ){ icPF = 0;}
    if( decay == 2 ){ icPF = 1.;}
  }  
}

  


if(icPF != 0){
// // //    Exchange of pion-like resonances

    Matrix AUX_pf;  
    
    for( int i=0; i<4; i++ ){

      if(ChPT == 1){

      AUX_pf = ( I * (-1.) * icPF * (gA*2.*MN/(sqrt(2.)*fpi)) * Fbpipi/(t-Mpi2) ) * Gamma5 ;

//       AUX_pf = ( I * (-1.) * icPF * (gA*2.*MN/(sqrt(2.)*fpi)) * wFCT/(t-Mpi2) ) * Gamma5 ;
     	  //kpi - k_t = kpi - (Q - kpi) = 2*kpi - Q 
	  PF_regge[i] = (2.*kpi[i]-Q[i]) * AUX_pf ;
	
      }else{
	
      AUX_pf = ( I * (-1.) * icPF * (sqrt(2.)*gPiNN) * Fbpipi/(t-Mpi2) ) * Gamma5 ;
      
	PF_regge[i] = (2.*kpi[i]-Q[i]) * AUX_pf ;
      }
      
    }
}



// // // // // // //   Exchange of nucleon-like resonances 
// // // // // // // // // // // // // // // // // // // //
      double icNP, icCNP; 

  if( process == 0 || process == 2 ){
    if( nucleon==1 ){ 
      if( decay==1 ){ // p(e,e'pi0)p
	icNP = sqrt(1./2.); 
	icCNP = sqrt(1./2.);
      } 
      if( decay==2 ){ // p(e,e'pi+)n 
	icNP = 1.; 
	icCNP = 1.;
      } 
    }
    if( nucleon == 2 ){ 
      if( decay==1 ){ // n(e,e'pi0)n
	icNP = -sqrt(1./2.);
	icCNP = -sqrt(1./2.); 
      } 
      if( decay==2 ){ // n(e,e'pi-)p 
	//test ?
	icNP = 1.; 
	icCNP = 1.; 
      } 
    }
  }
// // // 
  if( process == 1 ){
    if( Helicity == -1 ){
      if( nucleon==1 ){ // p(nu,mu-pi+)p 
	icNP = 0; 
	icCNP = 1.;
      } 
      if( nucleon == 2 ){ 
	if( decay==1 ){ // n(nu,mu-pi0)p 
	  icNP = sqrt(1./2.); 
	  icCNP = -sqrt(1./2.);   
	} 
	if( decay==2 ){ // n(nu,mu-pi+)n 
	  icNP = 1.;
	  icCNP = 0;   
	} 
      }
   }
   
    
//   HELICITY == 1 --> W^- induced 1-pion production 
// 	nucleon = 2 --> neutron initial state
// 		decay = 1 or 2 --> n + pi^-
// 	nucleon = 1 --> proton initial state
// 		decay = 1 --> n + pi^0
// 		decay = 2 --> p + pi^+

    if( Helicity == 1 ){
      if( nucleon==2 ){ //nub + n -->n + pi^-
	icNP = 0; 
	icCNP = 1.; 
      } 
      if( nucleon == 1 ){ 
	if( decay==1 ){ //nub + p --> n + pi^0
	  icNP = -sqrt(1./2.); 
	  icCNP = sqrt(1./2.); 
	} 
	if( decay==2 ){ //nub + p --> p + pi^-
	  icNP = 1.;
	  icCNP = 0;
	} 
      }
   }
  }


  Matrix AUX_np, AUX_cnp;
  
  Matrix AUX_np_V, AUX_cnp_V;
  Matrix AUX_np_A, AUX_cnp_A;
  
  Matrix QSlash_5;
    QSlash_5 = QSlash*Gamma5;
    
  if( icNP != 0 ){
    
    Matrix common;
    if( ChPT == 1 ){
      common = ( I*(-1.)*(icNP)*(gA/(sqrt(2.)*fpi))  / (s - MN2) ) * ( kpiSlash * Gamma5 * (sSlash + MN*Id) );
    }else{
      common = ( I*(icNP*sqrt(2.))*gPiNN  / (s - MN2) ) * ( Gamma5 * (sSlash + MN*Id) );
    }
    
    AUX_np_V = Fp * common;
    
    if( process == 1 || process == 2 ){
      AUX_np_A = GA * common;    
      for( int i=0; i<4; i++ ){
	NPv_regge[i] = AUX_np_V*Gamma[i];

	//ALEXIS: Can make choice to exclude NP from axial current in future:
	NPa_regge[i] = -1.*AUX_np_A*(Gamma_mu5[i] + (Q[i]/(Mpi2+Qsq))*QSlash_5); 
      }
    }
    
    if( process == 0 ){
      for( int i=0; i<4; i++ ){
	NPv_regge[i] = AUX_np_V*Gamma[i];
	NPa_regge[i] = 0.*Gamma5; 
      }
    }
    
  } //if( icNP != 0 )
  
  
  
  if( icCNP != 0 ){
    
    Matrix common_c;
    if( ChPT == 1 ){
      common_c = ( I*(-1.)*(icCNP)*(gA/(sqrt(2.)*fpi)) / (u - MN2) ) * ( (uSlash + MN*Id) * kpiSlash * Gamma5  );
    }else{
      common_c = ( I*(icCNP*sqrt(2.))*gPiNN / (u - MN2) ) * ( (uSlash + MN*Id) * Gamma5  );
    }
    
    AUX_cnp_V = Fp_c * common_c ;

    if( process == 1 || process == 2 ){
      AUX_cnp_A = GA_c * common_c ;
      for( int i=0; i<4; i++ ){
	CNPv_regge[i] = Gamma[i] * AUX_cnp_V ; 

	//ALEXIS: Can make choice to exclude CNP from axial current in future:
	CNPa_regge[i] =  -1.*(Gamma_mu5[i] + (Q[i]/(Mpi2+Qsq))*QSlash_5)*AUX_cnp_A; 
      }
    }
    if( process == 0 ){
      for( int i=0; i<4; i++ ){
	CNPv_regge[i] = Gamma[i] * AUX_cnp_V /*+ (Fp_c*I*(icCNP*sqrt(2.)/Qsq)*gPiNN)*(-Q[i]) * Gamma5*/; 
	CNPa_regge[i] = 0*Gamma5;
      }
    }
  }
  
  
  
if( ChPT == 1 && icPF != 0 ){ //Contact Term 
    
    complex<double> AUX_CT;  
    AUX_CT = ( I * (-1.) * icPF*(1./(sqrt(2.)*fpi)) );  
    
    double krho = 0.;
    double mrho2, Fwrhopi, Frho;
    complex<double> AUX_CT_A;
    Matrix AUX_PP;
    Matrix sigmamunu_Qt;
  
    if( process == 1 || process == 2){
      
      mrho2 = pow(775.8,2);
      double lambda_rho = 1260.;
      Fwrhopi = 1./( 1.+Qsq/pow(lambda_rho,2) );
      Frho = 1./(1.-t/mrho2)  * Fwrhopi ;
      
//       AUX_CT_A = (-1.)*AUX_CT*Frho*PRegge_NP;
      AUX_CT_A = (-1.)*AUX_CT*Frho*PRegge_rho;
      
      sigmamunu_Qt = QSlash*tSlash-tSlash*QSlash; 
// 	AUX_PP = ( PRegge_NP * AUX_CT * Frho / (-Qsq-Mpi2) ) *QSlash;
// 	AUX_PP = ( PRegge_rho * AUX_CT * Frho / (-Qsq-Mpi2) ) *QSlash;
      AUX_PP = ( PRegge_rho * AUX_CT * Frho / (-Qsq-Mpi2) ) * (QSlash + (-krho/(4.*MN)) * sigmamunu_Qt ) ;
    }  

    
  for( int i=0; i<4; i++ ){
    
      if( process == 0 ){
	if(nucleon == 1 && decay == 2){
	  CTv_regge[i] = (AUX_CT * gA * Fp) * Gamma_mu5[i];
	}
	if(nucleon == 2 && decay == 2){
	  CTv_regge[i] = (AUX_CT * gA * Fp_c) * Gamma_mu5[i];
	}
	if(decay == 1){ //CT does not enter for Pi0
	  CTv_regge[i] = 0.*Id;
	}
      }

      if( process == 2 ){ //respect to process==0, we include axial CT and PP
	
	if(decay == 2){
	  CTv_regge[i] = (AUX_CT * gA * wFCT) * Gamma_mu5[i]; 
	  CTa_regge[i] = AUX_CT_A* ( Gamma[i] + ( -krho/(4.*MN) )* (Gamma[i]*tSlash - tSlash*Gamma[i]) ) ;
	  PP_regge[i]  = Q[i]*AUX_PP;      
	}
	
	if(decay == 1){//CT does not enter for Pi0
	  CTv_regge[i] = 0.*Id;
	  CTa_regge[i] = 0.*Id;
	  PP_regge[i]  = 0.*Id;
	}
	
      }
    
    
      if( process == 1 ){
	  
	if( Helicity == -1 ){
	  if( nucleon==1 ){ // p(nu,mu-pi+)p 
    // 	icNP = 0; 
    // 	icCNP = 1.;
	    CTv_regge[i] = (AUX_CT * gA * Fp_c) * Gamma_mu5[i];
	    CTa_regge[i] =  AUX_CT_A* ( Gamma[i] + ( -krho/(4.*MN) )* (Gamma[i]*tSlash - tSlash*Gamma[i]) ) ;
	    PP_regge[i]  =  Q[i]*AUX_PP;      
	  } 
	  if( nucleon == 2 ){ 
	    if( decay==1 ){ // n(nu,mu-pi0)p 
    // 	  icNP = sqrt(1./2.); 
    // 	  icCNP = -sqrt(1./2.);
	      CTv_regge[i] = (AUX_CT * gA * 0.5*(Fp+Fp_c)) * Gamma_mu5[i];
	      CTa_regge[i] = AUX_CT_A* ( Gamma[i] + ( -krho/(4.*MN) )* (Gamma[i]*tSlash - tSlash*Gamma[i]) ) ;
	      PP_regge[i]  = Q[i]*AUX_PP;
	    } 
	    if( decay==2 ){ // n(nu,mu-pi+)n 
    // 	  icNP = 1.;
    // 	  icCNP = 0;   
	    CTv_regge[i] = (AUX_CT * gA * Fp) * Gamma_mu5[i];
	    CTa_regge[i] = AUX_CT_A* ( Gamma[i] + ( -krho/(4.*MN) )* (Gamma[i]*tSlash - tSlash*Gamma[i]) ) ;
	    PP_regge[i]  = Q[i]*AUX_PP;
	    } 
	  }
      }
	if( Helicity == 1 ){
	  if( nucleon==2 ){ 
    // 	icNP = 0; 
    // 	icCNP = 1.; 
	    CTv_regge[i] = (AUX_CT * gA * Fp_c) * Gamma_mu5[i] ;
	    CTa_regge[i] =  AUX_CT_A* ( Gamma[i] + ( -krho/(4.*MN) )* (Gamma[i]*tSlash - tSlash*Gamma[i]) );
	    PP_regge[i]  =  Q[i]*AUX_PP;
	  } 
	  if( nucleon == 1 ){ 
	    if( decay==1 ){ 
    // 	  icNP = -sqrt(1./2.); 
    // 	  icCNP = sqrt(1./2.); 
	      CTv_regge[i] = (AUX_CT * gA * 0.5*(Fp+Fp_c)) * Gamma_mu5[i];
	      CTa_regge[i] =  AUX_CT_A* ( Gamma[i] + ( -krho/(4.*MN) )* (Gamma[i]*tSlash - tSlash*Gamma[i]) );
	      PP_regge[i]  =  Q[i]*AUX_PP;
	    } 
	    if( decay==2 ){ 
    // 	  icNP = 1.;
    // 	  icCNP = 0;
	      CTv_regge[i] = (AUX_CT * gA * Fp) * Gamma_mu5[i];
	      CTa_regge[i] =  AUX_CT_A* ( Gamma[i] + ( -krho/(4.*MN) )* (Gamma[i]*tSlash - tSlash*Gamma[i]) );
	      PP_regge[i]  =  Q[i]*AUX_PP;
	    } 
	  }
	}
      } //if( process == 1 )
    
  }// for( int i=0; i<4; i++ )
  
} //if( ChPT == 1 && icPF != 0 ){ //Contact Term 
  
 
}


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


#endif











