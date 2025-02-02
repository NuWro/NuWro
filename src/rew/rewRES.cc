#include "../event1.h"
#include "../nucleus.h"
#include "../params.h"
#include "ff.h"
#include "rewparams.h"

#include "../dis/LeptonMass.h"
#include "../dis/alfa.h"
#include "../dis/charge.h"
#include "../dis/delta.h"
#include "../dis/dis_cr_sec.h"
#include "../dis/res_kinematics.h"
#include "../dis/res_xsec.h"
#include "../dis/resevent2.h"
#include "../dis/singlepion.h"

// extern double SPP[2][2][2][3][40];

extern "C" {
void shhpythiaitokay_(void);
void youcanspeaknowpythia_(void);
}

void SetupSPP(params &param) {
  if (true) {  //! CheckSPPSetup()){ -- No way to know in general
    std::cout << "[INFO]: Setting up singlepion tables..." << std::endl;
    shhpythiaitokay_();
    singlepion(param);
    youcanspeaknowpythia_();
    std::cout << "[INFO]: Set up singlepion tables!" << std::endl;
  }
}

void SetupSPP() {
  std::cout << "[WARN]: Setting up singlepion tables with default parameter "
               "set."
            << std::endl;
  params default_p;
  SetupSPP(default_p);
}

/*
 * find PDG codes for final state hadrons
 * returns false for non-SPP events
 */
bool get_pdg_spp(event &e, int &pion_pdg, int &nucleon_pdg) {
  // first look for pion and nucleon
  for (unsigned int i = 0; i < e.out.size(); ++i)
    if (e.out[i].pion())
      pion_pdg = e.out[i].pdg;
    else if (e.out[i].nucleon())
      nucleon_pdg = e.out[i].pdg;

  // SPP events have lepton, pion, and nucleon in the primary vertex
  if (e.out.size() != 3) return false;

  return pion_pdg != 0 and nucleon_pdg != 0;
}

// get the fraction of created pi
double get_pi_fraction(const int pion_pdg, res_xsec xsec) {
  switch (pion_pdg) {
    case PDG::pdg_piP:
      return (xsec.delta_pip + xsec.dis_pip) / (xsec.delta_total + xsec.dis_total);
    case PDG::pdg_pi:
      return (xsec.delta_pi0 + xsec.dis_pi0) / (xsec.delta_total + xsec.dis_total);
    case -PDG::pdg_piP:
      return (xsec.delta_pim + xsec.dis_pim) / (xsec.delta_total + xsec.dis_total);
    default:
      return 1;
  }
}

double calcRES(event &e, params &p, nucleus &t) {
  if (!e.flag.res) return 1;

  res_kinematics kin(e, t);
  kin.set_kinematics(e);

  res_xsec xsec(kin, e.flag.cc);

  int pion_pdg = 0, nucleon_pdg = 0;                    // placeholders for final hadrons PDG
  bool is_spp = get_pdg_spp(e, pion_pdg, nucleon_pdg);  // SPP above Pythia threshold

  double pi_scaling = 1.0;  // scale given pi contribution

  if (not kin.is_above_pythia_threshold() || xsec.is_no_dis()) {  // below Pythia threshold
    xsec.set_xsec_nopythia(kin, p);
    pi_scaling = get_pi_fraction(pion_pdg, xsec);
  } else if (is_spp)  // single pion production
    xsec.set_xsec(kin, p, pion_pdg, nucleon_pdg, e.in[0].t);

  // for more inelastic or kaon production there is no need to call set_xsec

  return e.res_angrew * xsec.get_total(kin.jacobian) * pi_scaling;
}
