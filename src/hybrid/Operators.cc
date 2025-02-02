#include "Operators.h"
#include "diracMatrices.h"
#include "Delta_free.h"
#include "D13.h"
#include "S11.h"
#include "P11.h"
#include "NON-Resonant.h"
#include "CT_PF_PP.h"
#include "Regge_CPT.h"


// // // // // // // // // Operators which do not depend on "r" // // // // // // // // // 
void GetOperator(Hadron_prime* kin, Reaction_parameters* par,  Matrix Operator[])
{
   	double W_up  = 2400; //Above W_up only Regge is used
	double W_low = 1300; //Below 1300 no Regge is added

	//Get low-W BG and Resonances    	
	Matrix BG[4], Resonances[4];
	if (kin->W < W_up){
            Get_Res_ChPT(kin, par, BG, Resonances);
	}


	//Get Reggeized Background
	Matrix ReChi[4];
	if (kin->W > W_low){
	    if( (-1.*kin->t) < kin->s  && -1.*kin->t < 10.*1e6) 
	    {
		    Get_ReChi(kin, par, ReChi);
	    }
	}

	//Mixing the BG and Regge
	
        double W0 = 1.5E3; // MeV, center of the transition
        double del_W0 = 0.1E3; // MeV, half width of the transition 

        double Phi_mix = Pi/2.*( 1.- 1./( 1.+exp( (kin->W-W0)/del_W0 ) ) );
        double cos2 = pow(cos(Phi_mix),2);
        double sin2 = pow(sin(Phi_mix),2);


        for( int i=0; i<4; i++ )
	{
	    for (int l = 0 ; l < 4 ; l++)
	    {
		for (int k = 0 ; k < 4 ; k++)
		{
		  Operator[i].M[l][k] = Resonances[i].M[l][k] + cos2*BG[i].M[l][k] + sin2*ReChi[i].M[l][k]; // Hybrid model
	    	}
	    }
        }

}


void Get_ReChi(Hadron_prime* kin, Reaction_parameters* par, Matrix ReChi[])
{

      //Set parameters:
      int process = par->process;
      int nucleon = par->nucleon;
      int decay = par->decay;
      int Helicity = par->Helicity;


       Matrix PF_regge[4], NPv_regge[4], NPa_regge[4], CNPv_regge[4], CNPa_regge[4], CTa_regge[4],
	 CTv_regge[4], PP_regge[4];

        Regge( process, nucleon, decay, Helicity, kin->Qsq, kin->ki, kin->kN, kin->Q, kin->kpi, kin->QSlash, kin->kpiSlash, kin->u, kin->s, kin->t, kin->sSlash, kin->uSlash, kin->tSlash,
	 PF_regge,  NPv_regge, NPa_regge, CNPv_regge, CNPa_regge, CTa_regge, CTv_regge, PP_regge );

	//The vector axial separation is not neccesary here so
	double multV=1.;
	double multA=1.;

        for( int i=0; i<4; i++ ){
	  for (int l = 0 ; l < 4 ; l++)
	  {
	    for (int k = 0 ; k < 4 ; k++)
	    {
          	ReChi[i].M[l][k] = multV*PF_regge[i].M[l][k] + multV*CTv_regge[i].M[l][k] + multA*CTa_regge[i].M[l][k]
				 + multV*NPv_regge[i].M[l][k] + multA*NPa_regge[i].M[l][k]
				 + multV*CNPv_regge[i].M[l][k]+ multA*CNPa_regge[i].M[l][k]
				 + multA*PP_regge[i].M[l][k];
	    }
	  }
	  
//             cout << "ReChi= " << ReChi[i] << endl;
        
        }

}


void Get_Res_ChPT(Hadron_prime* kin, Reaction_parameters* par, Matrix ChPT[], Matrix Resonances[])
{
      //Set parameters:
      int process = par->process;
      int nucleon = par->nucleon;
      int decay = par->decay;
      int Helicity = par->Helicity;
      int medmod = par->medmod; //OSMM Only for interaction with the nucleus, here medmod set to zero to make sure
      medmod = 0;

      //Some kinematic invariants:
      double s = kin->s;
      double t = kin->t;
      double u = kin->u;
      double Qsq = kin->Qsq;

      int cross;
      // ///////////////////////////////////
      // // Non-resonant contributions // // 
      /////////////////////////////////////

      // // NUCLEON POLE (non-resonant)
      cross = 0;
      Matrix Op_NP[4];
        NP_current( process, nucleon, decay, Helicity, cross, Qsq, s, kin->Q, kin->QSlash, kin->sSlash, kin->kpiSlash, Op_NP );


      // // CROSS NUCLEON POLE (Non-resonant)
      cross = 1;
      Matrix Op_NP_cross[4];
        CNP_current( process, nucleon, decay, Helicity, cross, Qsq, u, kin->Q, kin->QSlash, kin->uSlash, kin->kpiSlash, Op_NP_cross );

       Matrix Op_CT[4], Op_PF[4], Op_PP[4];
        CT_PF_PP( process, nucleon, decay, Helicity, kin->Q, kin->kpi, t, Qsq, kin->QSlash, Op_CT, Op_PF, Op_PP );

      // // // // // // // // // // // // // // // // // // // // // // // // // // // // //

      // // HIGHER RESONANCES       
      cross = 0;
      Matrix Op_delta[4];
      int Pascalutsa = 0; // Pascalutsa = 1 for Pascalutsa decay vertex, needs to be checked -> use 0
      DP_current( medmod, Pascalutsa, process, nucleon, decay, Helicity, cross, Qsq, s, u, kin->Q, kin->sMan, kin->ki, kin->kpi, kin->sSlash, Op_delta );
      double Deltaff;
      Delta_ff(s, u, Deltaff);

// // D13 POLE
      cross = 0;
      Matrix Op_D13[4];
      D13P_current( process, nucleon, decay, Helicity, cross, Qsq, s, kin->Q, kin->sMan, kin->ki, kin->kpi, kin->sSlash, Op_D13 );
      double D13ff;
      D13_ff(s, u, D13ff);

      // // // P11 POLE
      cross = 0;
      Matrix Op_P11[4];
      P11P_current( process, nucleon, decay, Helicity, cross, Qsq, s, kin->Q, kin->sMan, kin->kpiSlash, kin->QSlash, kin->sSlash, Op_P11 );
      double P11ff;
      P11_ff(s, u, P11ff);

      // // // // S11 POLE
      cross = 0;
      Matrix Op_S11[4];
      S11P_current( process, nucleon, decay, Helicity, cross, Qsq, s, kin->Q, kin->sMan, kin->kpiSlash, kin->QSlash, kin->sSlash, Op_S11 );
      double S11ff;
      S11_ff(s, u, S11ff);

      // // // // // // // // // // // // // // // // // // // // // // // // // // //

      // // HIGHER RESONANCES (u-channel)
      //       
      // // CROSS D13 POLE
      cross = 1;
      Matrix Op_D13_cross[4];
        CD13P_current( process, nucleon, decay, Helicity, cross, Qsq, u, kin->minusQ, kin->uMan, kin->kN, kin->kpi, kin->uSlash, Op_D13_cross );


      // // // // CROSS P11 POLE
      cross = 1;
      Matrix Op_P11_cross[4];
        CP11P_current( process, nucleon, decay, Helicity, cross, Qsq, u, kin->Q, kin->uMan, kin->kpiSlash, kin->QSlash,kin->uSlash,Op_P11_cross );

      // // // // CROSS S11 POLE
      cross = 1;
      Matrix Op_S11_cross[4];
        CS11P_current( process, nucleon, decay, Helicity, cross, Qsq, u, kin->Q, kin->uMan, kin->kpiSlash, kin->QSlash, kin->uSlash, Op_S11_cross );

// // // CROSS DELTA POLE        
       cross = 1;
      Matrix Op_delta_cross[4];
        CDP_current( medmod, Pascalutsa, process, nucleon, decay, Helicity, cross, Qsq, s, u, kin->minusQ, kin->uMan, kin->kN, kin->kpi, kin->uSlash, Op_delta_cross );

      // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
      //
      //Hadronic form factors are symmetric in s and u

        double CD13ff =D13ff;
        double CS11ff=S11ff;
        double CP11ff=P11ff;
        double CDeltaff=Deltaff;

        
      for( int i=0; i<4; i++ ){
          for (int l = 0 ; l < 4 ; l++)
          {
            for (int k = 0 ; k < 4 ; k++)
            {
		ChPT[i].M[l][k] = Op_NP[i].M[l][k] 
				+ Op_NP_cross[i].M[l][k] 
				+ Op_PP[i].M[l][k] 
				+ Op_CT[i].M[l][k] 
				+ Op_PF[i].M[l][k];  //Non-resonant

		//Add also cross channel resonances to 'BG'
		ChPT[i].M[l][k] += (CD13ff*Op_D13_cross[i].M[l][k] 
				+   CS11ff*Op_P11_cross[i].M[l][k]
				+   CP11ff*Op_S11_cross[i].M[l][k]
                                +   CDeltaff*Op_delta_cross[i].M[l][k]); 

		//Resonant contributions
 		Resonances[i].M[l][k] =   D13ff*Op_D13[i].M[l][k] 
                                 	+ P11ff*Op_P11[i].M[l][k] 
                                 	+ S11ff*Op_S11[i].M[l][k] 
                                  	+ Deltaff*Op_delta[i].M[l][k];


            } //k
          } //l
      } //i

    
    
}
