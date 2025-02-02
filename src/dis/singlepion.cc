#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
//#include "hist.h"
#include "pdg_name.h"
#include "LeptonMass.h"
#include "parameters.h"
#include "jednostki.h"
#include "masses.h"
//#include "bodek_old.h"
#include "grv94_bodek.h"
//#include "wiazka.h"
#include "dis_cr_sec.h"
#include "fragmentation.h"
#include "vect.h"
#include "charge.h"
//#include "lorentz.h"
#include <TMCParticle.h>
#include <TPythia6.h>
#include "singlepionhadr.h"
#include "params.h"
#include "singlepion.h"


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///// Computation of "single-pion functions" measuring fraction of SPP channels in overall DIS final states  /////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////  PYTHIA is run several times and the output information is put to tables: COUNTER, OVERALL SPP   ///////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////  The format is SPP[2][2][2][2][40] and the same for COUNTER and OVERALL ////////////////////////////////////////
/////  SPP[0]... neutrino
/////  SPP[1]... antineutrino
/////  SPP[ ][0]... CC reaction
/////  SPP[ ][1]... NC reaction
/////  SPP[ ][ ][0]... proton target
/////  SPP[ ][ ][1]... neutron target
/////  SPP[ ][ ][ ][0]... the channel with pi+ in the final state 
/////  SPP[ ][ ][ ][1]... the channel with pi0 in the final state 
/////  SPP[ ][ ][ ][2]... the channel with pi- in the final state 
/////  SPP[ ][ ][ ][ ][n].. invariant hadronic mass from 1210 to 1990 every 20 MeV
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////   SPP is then used to define actual spp function in the file alfa.cc
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double COUNTER[2][2][2][3][NSPPbins];
double OVERALL[2][2][2][NSPPbins];	//in OVERALL we do not distiguish channels
double SPP[2][2][2][3][NSPPbins];

double sppweight;

void
singlepion (params & p)		//produce SPP table
{
  double M2 = M12 * M12;
  int ile = p.spp_precision;


//Setting all the initial valus to be zero
  for (int j = 0; j < 2; j++)	//neutrino or antineutrino
    {
      for (int k = 0; k < 2; k++)	//CC or NC
	{
	  for (int l = 0; l < 2; l++)	//proton or neutron
	    {
	      for (int n = 0; n < 3; n++)	//channel choice
		{
		  for (int sa = 0; sa < NSPPbins; sa++)	//invariant hadronic mass
		    {
		      OVERALL[j][k][l][sa] = 0;
		      COUNTER[j][k][l][n][sa] = 0;
		    }
    }}}}

  double E = 20000;		//typical (anti-)neutrino energy; SPP should not depend on the choice of E
  double m;

  for (int j = 0; j < 2; j++)	//neutrino or antineutrino
    {
      for (int k = 0; k < 2; k++)	//CC or NC
	{
	  for (int l = 0; l < 2; l++)	//proton or neutron
	    {
	      for (int s = 0; s < NSPPbins; s++)
		{
		  double W = SPP_MIN + s * NSPPSize;
		  if (k == 0)
		    m = 105;
		  if (k == 1)
		    m = 0;
          double m2 = m * m;
		  double nulow =
		    ((M12 + E) * (W * W - M2 - m2) + 2 * M12 * E * E -
		     E * sqrt (kwad (W * W - M2 - m2 - 2 * M12 * E) -
			       4 * m2 * M12 * (M12 +
					       2 * E))) / 2 / M12 / (M12 +
								     2 * E) +
		    5;
		  double nuhigh =
		    ((M12 + E) * (W * W - M2 - m2) + 2 * M12 * E * E +
		     E * sqrt (kwad (W * W - M2 - m2 - 2 * M12 * E) -
			       4 * m2 * M12 * (M12 +
					       2 * E))) / 2 / M12 / (M12 +
								     2 * E) -
		    5;

		  {
		    for (double nu = nulow; nu < nuhigh; nu += (nuhigh - nulow) / ile)	//average over several values of the energy transfer
		      {

			int lepton;
			if (j == 0)
			  lepton = 14;
			else
			  lepton = -14;

			bool current;
			if (k == 0)
			  current = true;
			else
			  current = false;

			double m = lepton_mass (abs (lepton), current);	//mass of the given lepton (see pgd header file)

			int nukleon;
			if (l == 0)
			  nukleon = 2212;
			else
			  nukleon = 2112;

			sppweight =
			  cr_sec_dis (E, W, nu, lepton, nukleon, current);

			if (sppweight < 0.0)
			  continue;

			singlepionhadr (E, j, k, l, s, nu);

		      }
		  }
		}
	    }

	}
    }

//Setting SPP values
  for (int j = 0; j < 2; j++)	//neutrino or antineutrino
    {
      for (int k = 0; k < 2; k++)	//CC or NC
	{
	  for (int l = 0; l < 2; l++)	//proton or neutron
	    {
	      for (int n = 0; n < 3; n++)	//channel choice
		{
		  for (int s = 0; s < NSPPbins; s++)	//invariant hadronic mass
		    {
		      if (OVERALL[j][k][l][s] == 0)
			SPP[j][k][l][n][s] = 0;
		      else
			{
			  SPP[j][k][l][n][s] =
			    double (COUNTER[j][k][l][n][s] /
				    OVERALL[j][k][l][s]);
			}
		    }
		}
	    }
	}
    }
}

//koniec
