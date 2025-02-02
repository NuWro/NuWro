#include "res_xsec.h"
#include "../rew/rewparams.h"
#include "alfa.h"
#include "delta.h"
#include "dis2res.h"
#include "dis_cr_sec.h"
#include "resevent2.h"
#include "singlepion.h"
#include "resevent_hybrid.h"


extern "C" int pycomp_(const int *);

/*double pdd_red(double energy) {
  if (energy >= 1000)
    return 0.85;
  else if (energy > 750)
    return 0.65 + energy * 0.05 / 250.0;
  else  // if (en<=750)
    return 0.2 + energy * 0.2 / 250.0;
}*/

res_xsec::res_xsec(res_kinematics &kin, const bool cc)
    : dis_pip(0),
      dis_pi0(0),
      dis_pim(0),
      delta_pip(0),
      delta_pi0(0),
      delta_pim(0),
      dis_total(0),
      delta_total(0),
      is_cc(cc) {
  // determine indices for SPP table (see singlepion.cc)
  j = kin.neutrino.pdg < 0;
  k = not cc;
  l = kin.target.pdg != PDG::pdg_proton;

  // total electric charge of the pion-nucleon system
  final_charge = charge(kin.target.pdg) + (1 - k) * (1 - 2 * j);

  // the contribution to the cross section coming from DIS
  from_dis = max(0.0, cr_sec_dis(kin.neutrino.E(), kin.W, kin.q.t, kin.neutrino.pdg, kin.target.pdg, is_cc));
}

inline double res_xsec::get_dis_spp(const int pion_code, const res_kinematics &kin, const params &p) {
  return from_dis * SPP[j][k][l][pion_code][0] * betadis(j, k, l, pion_code, kin.W, p.bkgrscaling, p.res_dis_blending_start, p.res_dis_blending_end);
}

inline double res_xsec::get_delta_spp(const int pion_pdg, const int nucleon_pdg, const res_kinematics &kin,
                                      const params &p) {
  return cr_sec_delta(p.delta_FF_set, rew.pion_axial_mass.val, rew.pion_C5A.val, kin.neutrino.t, kin.W,
                      kin.q.t, kin.effective_mass, kin.neutrino.pdg, kin.target.pdg, nucleon_pdg, pion_pdg, is_cc) *
         alfadelta(j, k, l, pdg2spp(pion_pdg), kin.W, p.res_dis_blending_start, p.res_dis_blending_end);
}

int res_xsec::get_pion_pdg() {
  const double pip_fraction = (dis_pip + delta_pip) / (dis_total + delta_total);
  const double pi0_fraction = (dis_pi0 + delta_pi0) / (dis_total + delta_total);

  // randomly select final state pion
  double rand01 = frandom();

  if (pip_fraction > rand01)
    return PDG::pdg_piP;
  else if (pip_fraction + pi0_fraction > rand01)
    return PDG::pdg_pi;
  else
    return -PDG::pdg_piP;
}

/*
 * cross section calculated from SPP tables and Delta
 * it is used if invariant mass is smaller than 1210 (pythia threshold set in res_kinematics.cc)
 * or the dis contribution is determined to be 0 (from_dis)
 */
void res_xsec::set_xsec_nopythia(const res_kinematics &kin, const params &p) {
  // contributions from DIS
  dis_pip = get_dis_spp(pip, kin, p);
  dis_pi0 = get_dis_spp(pi0, kin, p);
  dis_pim = get_dis_spp(pim, kin, p);

  dis_total = dis_pip + dis_pi0 + dis_pim;

  // contributions from Delta
  switch (final_charge) {
    case 2:  // pi+ + proton
      delta_pip = get_delta_spp(PDG::pdg_piP, PDG::pdg_proton, kin, p);
      break;
    case 1:  // pi+ + neutron or pi0 + proton
      delta_pip = get_delta_spp(PDG::pdg_piP, PDG::pdg_neutron, kin, p);
      delta_pi0 = get_delta_spp(PDG::pdg_pi, PDG::pdg_proton, kin, p);
      break;
    case 0:  // pi0 + neutron or pi- + proton
      delta_pi0 = get_delta_spp(PDG::pdg_pi, PDG::pdg_neutron, kin, p);
      delta_pim = get_delta_spp(-PDG::pdg_piP, PDG::pdg_proton, kin, p);
      break;
    case -1:  // pi- + neutron
      delta_pim = get_delta_spp(-PDG::pdg_piP, PDG::pdg_neutron, kin, p);
      break;
    default:
      cerr << "[WARNING]: charge out of rangen\n";
  };

  delta_total = delta_pip + delta_pi0 + delta_pim;
}

// cross section calculated based on generated final states from Pythia
void res_xsec::set_xsec(res_kinematics &kin, const params &p, const int pion_pdg, const int nucleon_pdg,
                        const double neutrino_energy) {
  // PDG to SPP code
  const int t = pdg2spp(pion_pdg);

  // dis contribution to single pion production
  dis_total = from_dis * alfadis(j, k, l, t, kin.W, p.res_dis_blending_start, p.res_dis_blending_end);

  // delta contribution to single pion production
  delta_total = cr_sec_delta(p.delta_FF_set, rew.pion_axial_mass.val, rew.pion_C5A.val, kin.neutrino.E(), kin.W,
                             kin.q.t, kin.effective_mass, kin.neutrino.pdg, kin.target.pdg, nucleon_pdg, pion_pdg, is_cc) /
                SPPF(j, k, l, t, kin.W) * alfadelta(j, k, l, t, kin.W, p.res_dis_blending_start, p.res_dis_blending_end);

  // reduce cross section by removing the contribution from pionless delta decay
  // more details in: J. Żmuda and J.T. Sobczyk, Phys. Rev. C 87, 065503 (2013)
  // if ((p.nucleus_p + p.nucleus_n) > 7) delta_total *= pdd_red(neutrino_energy);
}

void res_xsec_hybrid::set_xsec_nopythia(const res_kinematics &kin,
                                        const params &p, const double xsec_pip,
                                        const double xsec_pi0,
                                        const double xsec_pim) {
  // contributions from DIS
  // dis_pip = get_dis_spp(pip, kin, p);
  // dis_pi0 = get_dis_spp(pi0, kin, p);
  // dis_pim = get_dis_spp(pim, kin, p);

  dis_total = 0.;
  // contributions from Delta
  // auto pion_momentum = get_pion_momentum(kin.W);

  delta_pip = xsec_pip * alfadelta(0, 1, 0, pdg2spp(211), kin.W, p.res_dis_blending_start, p.res_dis_blending_end);
  delta_pi0 = xsec_pi0 * alfadelta(0, 1, 0, pdg2spp(111), kin.W, p.res_dis_blending_start, p.res_dis_blending_end);
  delta_pim = xsec_pim * alfadelta(0, 1, 0, pdg2spp(-211), kin.W, p.res_dis_blending_start, p.res_dis_blending_end);

  // delta_total = delta_pip + delta_pi0 + delta_pim;
  int final_charge = charge(kin.target.pdg) + (1 - k) * (1 - 2 * j);
  switch (final_charge) {
  case -1:
      delta_total = delta_pim;
      break;
  case 0:
      delta_total = delta_pi0 + delta_pim;
      break;
  case 1:
      delta_total = delta_pip + delta_pi0;
      break;
  case 2:
      delta_total = delta_pip;
      break;
  };
  // return delta_total;
}

void res_xsec_hybrid::set_xsec(res_kinematics &kin, const params &p,
                               const int pion_pdg, const int nucleon_pdg,
                               const double neutrino_energy,
                               const double hybrid_channel_xsec) {
  // PDG to SPP code
  const int t = pdg2spp(pion_pdg);

  // dis contribution to single pion production
  dis_total = from_dis * alfadis(0, 1, 0, t, kin.W, p.res_dis_blending_start, p.res_dis_blending_end);

  delta_total = hybrid_channel_xsec / SPPF(j, k, l, t, kin.W) *
                alfadelta(0, 1, 0, t, kin.W, p.res_dis_blending_start, p.res_dis_blending_end);
  // reduce cross section by removing the contribution from pionless delta decay
  // more details in: J. Żmuda and J.T. Sobczyk, Phys. Rev. C 87, 065503 (2013)
  // if ((p.nucleus_p + p.nucleus_n) > 7)
  //   delta_total *= pdd_red(neutrino_energy);
}

double res_xsec_hybrid::get_dis_spp(const int pion_code,
                                    const res_kinematics &kin,
                                    const params &p) {
  return from_dis * SPP[j][k][l][pion_code][0] *
         alfadis(0, 1, 0, pion_code, kin.W, p.res_dis_blending_start, p.res_dis_blending_end);
}
