#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "params.h"
//#include "hist.h"
#include "pdg_name.h"
#include "parameters.h"
#include "LeptonMass.h"
#include "jednostki.h"
//#include "bodek_old.h"
#include "grv94_bodek.h"
//#include "wiazka.h"
#include "dis_cr_sec.h"
//#include "delta.h"
#include "fragmentation.h"
//#include "vect.h"
//#include "charge.h"
//#include "lorentz.h"
//#include "dis2res.h"
#include "event1.h"
//#include "sobek.h"
#include <TMCParticle.h>
#include <TPythia6.h>
#include "singlepionhadr.h"
#include "singlepion.h"

TPythia6 *pythia3 = new TPythia6 ();
extern "C" int pycomp_ (const int *);

// extern double COUNTER[2][2][2][3][40];
// extern double OVERALL[2][2][2][40];	//in OVERALL we do not distiguish channels

extern double sppweight;

void
singlepionhadr (double E, int j, int k, int l, int s, double nu)
{
////////////////////////////////////////////////////////////////////////////
//      Setting Pythia parameters
//      Done by Jaroslaw Nowak
//////////////////////////////////////////////

//stabilne pi0
  pythia3->SetMDCY (pycomp_ (&pizero), 1, 0);

  pythia3->SetMSTU (20, 1);	//advirsory warning for unphysical flavour switch off
  pythia3->SetMSTU (23, 1);	//It sets counter of errors at 0
  pythia3->SetMSTU (26, 0);	//no warnings printed 

  // PARJ(32)(D=1GeV) is, with quark masses added, used to define the minimum allowable enrgy of a colour singlet parton system
  pythia3->SetPARJ (33, 0.1);

  // PARJ(33)-PARJ(34)(D=0.8GeV, 1.5GeV) are, with quark masses added, used to define the remaining energy below which 
  //the fragmentation of a parton system is stopped and two final hadrons formed.
  pythia3->SetPARJ (34, 0.5);
  pythia3->SetPARJ (35, 1.0);

  //PARJ(36) (D=2.0GeV) represents the dependence of the mass of final quark pair for defining the stopping point of the 
  //fragmentation. Strongly corlated with PARJ(33-35)

  pythia3->SetPARJ (37, 1.);	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

  //MSTJ(17) (D=2) number of attemps made to find two hadrons that have a combined mass below the cluster mass and thus allow 
  // a cluster to decay rather than collaps
  pythia3->SetMSTJ (18, 3);	//do not change

/////////////////////////////////////////////////
//              End of setting Pythia parameters
////////////////////////////////////////////////

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

  double W = SPP_MIN + NSPPSize * s;

  int nParticle = 0;
  int nCharged = 0;
  int NPar = 0;
  Pyjets_t *pythiaParticle;	//deklaracja event recordu
  double W1 = W / GeV;		//W1 w GeV-ach potrzebne do Pythii

  while (NPar < 5)
    {
      hadronization (E, W, nu, m, lepton, nukleon, current);
      pythiaParticle = pythia3->GetPyjets ();
      NPar = pythia3->GetN ();
    }

  OVERALL[j][k][l][s] += sppweight;;	//overall counter

  if (pythia3->GetN () == 5)	//this means that the reaction was SPP
    {
      if (lepton > 0 && current == true && nukleon == neutron)
	{
	  if (pythiaParticle->K[1][3] == 211
	      || pythiaParticle->K[1][4] == 211)
	    {
	      COUNTER[j][k][l][0][s] += sppweight;
	    }			// nu + n -> mu + n + (pi+)
	  if (pythiaParticle->K[1][3] == 111
	      || pythiaParticle->K[1][4] == 111)
	    {
	      COUNTER[j][k][l][1][s] += sppweight;
	    }			// nu + n -> mu + p + (pi0)
	  if (pythiaParticle->K[1][3] == -211
	      || pythiaParticle->K[1][4] == -211)
	    {
	      cout << "error1" << endl;
	    }
	}
      if (lepton > 0 && current == true && nukleon == proton)
	{
	  if (pythiaParticle->K[1][3] == 211
	      || pythiaParticle->K[1][4] == 211)
	    {
	      COUNTER[j][k][l][0][s] += sppweight;
	    }			// nu + p -> mu + p + (pi+)
	  if (pythiaParticle->K[1][3] == 111
	      || pythiaParticle->K[1][4] == 111)
	    {
	      cout << "error2" << endl;
	    }
	  if (pythiaParticle->K[1][3] == -211
	      || pythiaParticle->K[1][4] == -211)
	    {
	      cout << "error3" << endl;
	    }
	}
      if (lepton > 0 && current == false && nukleon == neutron)
	{
	  if (pythiaParticle->K[1][3] == -211
	      || pythiaParticle->K[1][4] == -211)
	    {
	      COUNTER[j][k][l][2][s] += sppweight;
	    }			// nu + n -> nu + p + (pi-)
	  if (pythiaParticle->K[1][3] == 111
	      || pythiaParticle->K[1][4] == 111)
	    {
	      COUNTER[j][k][l][1][s] += sppweight;
	    }			// nu + n -> nu + n + (pi0)
	  if (pythiaParticle->K[1][3] == 211
	      || pythiaParticle->K[1][4] == 211)
	    {
	      cout << "error4" << endl;
	    }
	}
      if (lepton > 0 && current == false && nukleon == proton)
	{
	  if (pythiaParticle->K[1][3] == 211
	      || pythiaParticle->K[1][4] == 211)
	    {
	      COUNTER[j][k][l][0][s] += sppweight;
	    }			// nu + p -> nu + n + (pi+)
	  if (pythiaParticle->K[1][3] == 111
	      || pythiaParticle->K[1][4] == 111)
	    {
	      COUNTER[j][k][l][1][s] += sppweight;
	    }			// nu + p -> nu + p + (pi0)
	  if (pythiaParticle->K[1][3] == -211
	      || pythiaParticle->K[1][4] == -211)
	    {
	      cout << "error5" << endl;
	    }
	}
      if (lepton < 0 && current == true && nukleon == proton)
	{
	  if (pythiaParticle->K[1][3] == -211
	      || pythiaParticle->K[1][4] == -211)
	    {
	      COUNTER[j][k][l][2][s] += sppweight;
	    }			// anu + p -> amu + p + (pi-)
	  if (pythiaParticle->K[1][3] == 111
	      || pythiaParticle->K[1][4] == 111)
	    {
	      COUNTER[j][k][l][1][s] += sppweight;
	    }			// anu + p -> amu + n + (pi0)
	  if (pythiaParticle->K[1][3] == 211
	      || pythiaParticle->K[1][4] == 211)
	    {
	      cout << "error6" << endl;
	    }
	}
      if (lepton < 0 && current == true && nukleon == neutron)
	{
	  if (pythiaParticle->K[1][3] == -211
	      || pythiaParticle->K[1][4] == -211)
	    {
	      COUNTER[j][k][l][2][s] += sppweight;
	    }			// anu + n -> amu + n + (pi-)
	  if (pythiaParticle->K[1][3] == 111
	      || pythiaParticle->K[1][4] == 111)
	    {
	      cout << "error7" << endl;
	    }
	  if (pythiaParticle->K[1][3] == 211
	      || pythiaParticle->K[1][4] == 211)
	    {
	      cout << "error8" << endl;
	    }
	}
      if (lepton < 0 && current == false && nukleon == neutron)
	{
	  if (pythiaParticle->K[1][3] == -211
	      || pythiaParticle->K[1][4] == -211)
	    {
	      COUNTER[j][k][l][2][s] += sppweight;
	    }			// anu + n -> anu + p + (pi-)
	  if (pythiaParticle->K[1][3] == 111
	      || pythiaParticle->K[1][4] == 111)
	    {
	      COUNTER[j][k][l][1][s] += sppweight;
	    }			// anu + n -> anu + n + (pi0)
	  if (pythiaParticle->K[1][3] == 211
	      || pythiaParticle->K[1][4] == 211)
	    {
	      cout << "error9" << endl;
	    }
	}
      if (lepton < 0 && current == false && nukleon == proton)
	{
	  if (pythiaParticle->K[1][3] == 211
	      || pythiaParticle->K[1][4] == 211)
	    {
	      COUNTER[j][k][l][0][s] += sppweight;
	    }			// anu + p -> anu + n + (pi+)
	  if (pythiaParticle->K[1][3] == 111
	      || pythiaParticle->K[1][4] == 111)
	    {
	      COUNTER[j][k][l][1][s] += sppweight;
	    }			// anu + p -> anu + p + (pi0)
	  if (pythiaParticle->K[1][3] == -211
	      || pythiaParticle->K[1][4] == -211)
	    {
	      cout << "error10" << endl;
	    }

	}

    }

}

//koniec
