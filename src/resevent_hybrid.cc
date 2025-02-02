/*
This function calculates RES events from the Hybrid model
*/

#include "resevent_hybrid.h"
#include "TDecompLU.h"
#include "dis/LeptonMass.h"
#include "dis/res_kinematics.h"
#include "dis/res_xsec.h"
#include "dis/resevent2.h"
#include "event1.h"
#include "hybrid/hybrid_gateway.h"
#include "dis/fragmentation.h"
#include "mmapio.h"
#include "nucleus.h"
#include "params.h"
#include "particle.h"
#include "vect.h"
#include <array>
#include <cmath>
#include <string>
#include <signal.h>

class hybrid_grid_mmapio {
private:
  std::array<mmapio<double>, 6> tabs;
  hybrid_grid_mmapio(std::string basepath)
      : tabs(
            {mmapio<double>{basepath + "hybrid_grid_dQ2dWdcth_11.dat", false},
             mmapio<double>{basepath + "hybrid_grid_dQ2dWdcth_22.dat", false},
             mmapio<double>{basepath + "hybrid_grid_dQ2dWdcth_21.dat", false},
             mmapio<double>{basepath + "hybrid_grid_dQ2dWdcth_2_11.dat", false},
             mmapio<double>{basepath + "hybrid_grid_dQ2dWdcth_2_22.dat", false},
             mmapio<double>{basepath + "hybrid_grid_dQ2dWdcth_2_21.dat",
                            false}}) {
    std::cout << "loading hybrid grid from " << basepath << std::endl;
  }
  hybrid_grid_mmapio(std::string basepath, int)
      : tabs(
            {mmapio<double>{basepath + "hybrid_grid_dQ2dW_11.dat", false},
             mmapio<double>{basepath + "hybrid_grid_dQ2dW_22.dat", false},
             mmapio<double>{basepath + "hybrid_grid_dQ2dW_21.dat", false},
             mmapio<double>{basepath + "hybrid_grid_dQ2dW_2_11.dat", false},
             mmapio<double>{basepath + "hybrid_grid_dQ2dW_2_22.dat", false},
             mmapio<double>{basepath + "hybrid_grid_dQ2dW_2_21.dat",
                            false}}) {
    std::cout << "loading 2d hybrid grid from " << basepath << std::endl;
  }
public:
  static hybrid_grid_mmapio &get(char const *basepath = "") { 
    static hybrid_grid_mmapio instance_dQ2dWdcth{std::string(basepath)};
    return instance_dQ2dWdcth;
  }
  static hybrid_grid_mmapio &get2d(char const *basepath = "") { 
    static hybrid_grid_mmapio instance_dQ2dW{std::string(basepath), 1};
    return instance_dQ2dW;
  }
  const mmapio<double>& operator[](int i) const { return tabs[i]; }
};  // an adapter to mmapio for the hybrid grid
    // can be used as in place replacement for the old hybrid_grid array


constexpr double L0() { return 1.; }
constexpr double L1(double x) { return x; }
constexpr double L2(double x) { return (3 * x * x - 1) / 2; }

double poly_interop_1d(std::array<double , 3>& points, double x){
  double a0 = (points[0] + 4 * points[1] + points[2]) / 6;
  double a1 = (points[2] - points[0]) / 2;
  double a2 = (points[0] - 2 * points[1] + points[2]) / 3;
  return a0 * L0() + a1 * L1(x) + a2 * L2(x);
}

double poly_interop_2d(std::array<std::array<double , 3>,3>& points_2d, double x, double y){
  std::array<double , 3> points;
  for (int i = 0; i < 3; i++) {
    points[i] = poly_interop_1d(points_2d[i], y);
  }
  return poly_interop_1d(points, x);
}

double poly_interop_3d(std::array<std::array<std::array<double , 3>,3>,3>& points_3d, double x, double y, double z){
  std::array<double , 3> points;
  for (int i = 0; i < 3; i++) {
    points[i] = poly_interop_2d(points_3d[i], y, z);
  }
  return poly_interop_1d(points, x);
}

double poly_interop_4d(std::array<std::array<std::array<std::array<double , 3>,3>,3>,3>& points_4d, double x, double y, double z, double t){
  std::array<double , 3> points;
  for (int i = 0; i < 3; i++) {
    points[i] = poly_interop_3d(points_4d[i], y, z, t);
  }
  return poly_interop_1d(points, x);
}

size_t index_calculator(size_t cth_index, size_t W_index, size_t Q_index,
                        size_t Q_bins) { // use different Q_bins value to
                                         // identify different part of Q table
  return cth_index * Wbin_hybrid * Q_bins + W_index * Q_bins + Q_index;
}

size_t index_calculator_2d(size_t W_index, size_t Q_index,
                        size_t Q_bins) { // use different Q_bins value to
                                         // identify different part of Q table
  return W_index * Q_bins + Q_index;
}

void resevent_hybrid(params &p, event &e, nucleus& t, bool cc) // free nucleon only!
{
  e.weight = 0;           // if kinematically forbidden

  res_kinematics kin(e, t);  // kinematics variables

  // check threshold for pion production (otherwise left e.weight = 0)
  if (not kin.is_above_threshold()) return;

  // generate random kinematics (return false in the case of impossible kinematics)
  if (not kin.generate_kinematics(p.res_dis_cut)) return;
  //if (not kin.generate_kinematics(100000, 1310)) return;

  // save final lepton (kin.lepton is in target rest frame so boost it first)
  particle final_lepton = kin.lepton;
  final_lepton.boost(kin.target.v());
  final_lepton.pdg = kin.neutrino.pdg + cc * (1 - 2.0 * (kin.neutrino.pdg > 0));

  e.out.push_back(final_lepton);

  // final state particles
  particle final_pion, final_nucleon;
  particle final_pion2, final_nucleon2; // needed for choosing a decay channel

  // selection of the final state

  // specify the total electric charge of the pion-nucleon system
  int j = kin.neutrino.pdg < 0;
  int k = not cc;
  int l = kin.target.pdg != PDG::pdg_proton;
  int final_charge = charge(kin.target.pdg) + (1 - k) * (1 - 2 * j);

  double xsec_pip=0, xsec_pi0=0, xsec_pim=0, xsec_inclusive=0; // cross sections

  // choose a random direction in CMS
  vec kierunek = rand_dir();
  // or fix the direction in the Adler frame (tests)
  //vec kierunek = -hybrid_dir_from_adler(0, frandom()*2*Pi-Pi, kin.neutrino, kin.lepton);

  // specify the params needed for ABCDE (note strange order!)
  int params[4];
  params[0] = 1;                                   // only CC for now
  params[2] = (1 - 2.0 * (kin.neutrino.pdg > 0));  // helicity
  // params[3] is the target nucleon, params[1] is the decay channel

  // cross section function
  double (*hybrid_xsec)(res_kinematics*, int*, vect, double);
  switch(p.res_hybrid_sampling)
  {
    case 1: hybrid_xsec = hybrid_dsdQ2dW_tab; break;
    case 2: hybrid_xsec = hybrid_dsdQ2dWdcth_tab; break;
    case 3: hybrid_xsec = hybrid_dsdQ2dWdcth; break;
    case 4: hybrid_xsec = hybrid_dsdQ2dWdOm; break;
    default:hybrid_xsec = hybrid_dsdQ2dWdOm; break;
  }
  if(p.res_hybrid_sampling == 1){
    e.flag.need_resample_dir = true;
    e.flag.need_resample_phi = true;
  }
  // if tabularized grids are going to be used, load them here
  if (p.res_hybrid_sampling != 4){
    hybrid_grid_mmapio::get(p.table_path.c_str());
    hybrid_grid_mmapio::get2d(p.table_path.c_str());
  }
  // switch from tabs to dsdQ2dWdcth above the limits
  if (p.res_hybrid_sampling < 3 &&
      (kin.W > Wmax_hybrid || -kin.q * kin.q > Q2max_hybrid || -kin.q * kin.q < Q2min_hybrid)) {
    hybrid_xsec = hybrid_dsdQ2dWdcth;
    e.flag.need_resample_dir = false;
    e.flag.need_resample_phi = true;
  }

  // pion momentum if masses were averaged
  double pion_momentum;
  bool ghent_calc_done = false;
  auto do_ghent_xsec_calc = [&]() {
    if(ghent_calc_done){
      return true;
    }
    ghent_calc_done = true;
    switch (final_charge) { // calculate the cross section with only the CMS
                            // variables
    case 2:                 // pi+ + proton (nu_11)
    {
      params[3] = 1;
      params[1] = 1;
      final_pion.set_pdg_and_mass(PDG::pdg_piP);
      final_nucleon.set_pdg_and_mass(PDG::pdg_proton);
      pion_momentum = kin1part(kin.W, final_nucleon.pdg, final_pion.pdg,
                               final_nucleon, final_pion, kierunek);
      xsec_pip = hybrid_xsec(&kin, params, final_pion, pion_momentum);
    }
      xsec_inclusive = xsec_pip;
      if (not(xsec_inclusive > 0))
        return false;
      break;
    case 1: // pi+ + neutron (nu_22) or pi0 + proton (nu_21)
    {
      params[3] = 2;
      params[1] = 2;
      final_pion.set_pdg_and_mass(PDG::pdg_piP);
      final_nucleon.set_pdg_and_mass(PDG::pdg_neutron);
      pion_momentum = kin1part(kin.W, final_nucleon.pdg, final_pion.pdg,
                               final_nucleon, final_pion, kierunek);
      xsec_pip = hybrid_xsec(&kin, params, final_pion, pion_momentum);
    }
      {
        params[3] = 2;
        params[1] = 1;
        final_pion2.set_pdg_and_mass(PDG::pdg_pi);
        final_nucleon2.set_pdg_and_mass(PDG::pdg_proton);
        pion_momentum = kin1part(kin.W, final_nucleon2.pdg, final_pion2.pdg,
                                 final_nucleon2, final_pion2, kierunek);
        xsec_pi0 = hybrid_xsec(&kin, params, final_pion2, pion_momentum);
      }
      xsec_inclusive = xsec_pip + xsec_pi0;
      if (not(xsec_inclusive > 0 && xsec_pip > 0 && xsec_pi0 > 0))
        return false;
      break;
    case 0: // pi0 + neutron (anu_11) or pi- + proton (anu_12)
    {
      params[3] = 1;
      params[1] = 1;
      final_pion.set_pdg_and_mass(PDG::pdg_pi);
      final_nucleon.set_pdg_and_mass(PDG::pdg_neutron);
      pion_momentum = kin1part(kin.W, final_nucleon.pdg, final_pion.pdg,
                               final_nucleon, final_pion, kierunek);
      xsec_pi0 = hybrid_xsec(&kin, params, final_pion, pion_momentum);
    }
      {
        params[3] = 1;
        params[1] = 2;
        final_pion2.set_pdg_and_mass(-PDG::pdg_piP);
        final_nucleon2.set_pdg_and_mass(PDG::pdg_proton);
        pion_momentum = kin1part(kin.W, final_nucleon2.pdg, final_pion2.pdg,
                                 final_nucleon2, final_pion2, kierunek);
        xsec_pim = hybrid_xsec(&kin, params, final_pion2, pion_momentum);
      }
      xsec_inclusive = xsec_pi0 + xsec_pim;
      if (not(xsec_inclusive > 0 && xsec_pi0 > 0 && xsec_pim > 0))
        return false;
      break;
    case -1: // pi- + neutron (anu_21)
    {
      params[3] = 2;
      params[1] = 2;
      final_pion.set_pdg_and_mass(-PDG::pdg_piP);
      final_nucleon.set_pdg_and_mass(PDG::pdg_neutron);
      pion_momentum = kin1part(kin.W, final_nucleon.pdg, final_pion.pdg,
                               final_nucleon, final_pion, kierunek);
      xsec_pim = hybrid_xsec(&kin, params, final_pion, pion_momentum);
    }
      xsec_inclusive = xsec_pim;
      if (not(xsec_inclusive > 0 && xsec_pim > 0))
        return false;
      break;
    default:
      cerr << "[WARNING]: Reaction charge out of range\n";
    };
    return true;
  };

  e.flag.res_delta = false;
  const double factor =
      (G * G * cos2thetac / 2) /
      (kin.neutrino.E() * kin.neutrino.E() / 2 / kin.effective_mass / 1.) / cm2 * 1e38;

  // 0: unset, random selection for final charge = 1, 0 channel
  // 211/111/-211: force final pion to be in certain channel
  auto gen_final_particles_hybrid = [&](int pdg_pion = 0) {
    switch (final_charge) { // calculate the cross section with only the CMS
                            // variables
    case 2:                 // pi+ + proton (nu_11)
      break;
    case 1: // pi+ + neutron (nu_22) or pi0 + proton (nu_21)
      if (pdg_pion == 211){
        params[3] = 2;
        params[1] = 2;
        break;
      }
      if (pdg_pion == 111) {
        final_pion = final_pion2;
        final_nucleon = final_nucleon2;
        params[3] = 2;
        params[1] = 1;
        break;
      }
      if ((xsec_pip / xsec_inclusive <
          frandom())) // random selection, switch to "2"
      {
        final_pion = final_pion2;
        final_nucleon = final_nucleon2;
        params[3] = 2;
        params[1] = 1;
      } else // make sure the params are okay for "1"
      {
        params[3] = 2;
        params[1] = 2;
      }
      break;
    case 0: // pi0 + neutron (anu_11) or pi- + proton (anu_12)
      if (pdg_pion == 111) {
        params[3] = 1;
        params[1] = 1;
        break;
      }
      if (pdg_pion == -211) {
        final_pion = final_pion2;
        final_nucleon = final_nucleon2;
        params[3] = 1;
        params[1] = 2;
        break;
      }
      if (xsec_pi0 / xsec_inclusive <
          frandom()) // random selection, switch to "2"
      {
        final_pion = final_pion2;
        final_nucleon = final_nucleon2;
        params[3] = 1;
        params[1] = 2;
      } else // make sure the params are okay for "1"
      {
        params[3] = 1;
        params[1] = 1;
      }
      break;
    case -1: // pi- + neutron (anu_21)
      break;
    default:
      cerr << "[WARNING]: Reaction charge out of range\n";
    };
    e.flag.res_delta = true;
    final_nucleon.p4() = final_nucleon.boost(kin.hadron_speed);
    final_nucleon.p4() = final_nucleon.boost(kin.target.v());

    final_pion.p4() = final_pion.boost(kin.hadron_speed);
    final_pion.p4() = final_pion.boost(kin.target.v());

    // save final state hadrons
    // warning: the order is essential
    e.out.push_back(final_pion);
    e.out.push_back(final_nucleon);
  };
  res_xsec_hybrid xsec(kin, cc);
  if ((!kin.is_above_pythia_threshold()) || xsec.is_no_dis() ||
      (!p.dyn_dis_cc)) {
    // xsec.set_xsec_nopythia(kin, p, xsec_pip * factor ,
    //                        xsec_pi0 * factor , xsec_pim * factor );
    // e.weight = xsec.get_total(kin.jacobian);
    if(!do_ghent_xsec_calc()){
      return;
    }
    e.weight = xsec_inclusive * (factor* 1e-38) * kin.jacobian ;
    gen_final_particles_hybrid();
  } else // W above pythia threshold and fromdis > 0
  {
    // the algorithm starts from the production of PYTHIA event
    TPythia6 *pythia71 = get_pythia();

    int nof_particles = 0;      // number of particles in the final state
    Pyjets_t *pythia_particles; // pythia particles placeholder

    // force at least 5 particles in the final state
    // TODO: including initial particles?
    while (nof_particles < 5) {
      hadronization(kin.neutrino.E(), kin.W, kin.q.t, kin.lepton_mass,
                    kin.neutrino.pdg, kin.target.pdg, cc);
      pythia_particles = pythia71->GetPyjets();
      nof_particles = pythia71->GetN();
    }

    /*
    There are three different possible outcomes:
      a) spp event nof_particles = 5
      b) more inelastic event nof_particles > 5, typically 7
      c) single kaon production; this causes technical complications because
    nof_particles = 5 also in this case
    */

    // single pion production -> 5 particles including a pion
    if (nof_particles == 5 && (PDG::pion(pythia_particles->K[1][3]) ||
                                PDG::pion(pythia_particles->K[1][4]))) {
      if (!do_ghent_xsec_calc()) {
        return;
      }
      // K[1][3] or K[1][4] is a pion, and only 1 pion
      const int pion_pdg = PDG::pion(pythia_particles->K[1][3])
                               ? pythia_particles->K[1][3]
                               : pythia_particles->K[1][4];
      // if(pion_pdg)
      // PDG to SPP code
      const int t = pdg2spp(pion_pdg);
      const int nucleon_pdg =
          xsec.final_charge + t == 1 ? PDG::pdg_neutron : PDG::pdg_proton;

      // initialize xsec struct according to Pythia result
      double xsec_channel{};
      if (pion_pdg == 211){
        xsec_channel = xsec_pip*factor;
      }else if (pion_pdg == -211){
        xsec_channel = xsec_pim*factor;
      }else{
        xsec_channel = xsec_pi0*factor;
      }
      xsec.set_xsec(kin, p, pion_pdg, nucleon_pdg, e.in[0].t,xsec_channel);
      // xsec.set_xsec(kin, p, pion_pdg, nucleon_pdg, e.in[0].t,xsec_inclusive*factor);
      // xsec.set_delta_total(xsec_inclusive*factor);

      e.weight = xsec.get_total(kin.jacobian); // save total cross sections

      // randomly decide if SPP comes from Delta or DIS
      if (xsec.get_dis_fraction() > frandom()) { // SPP from DIS
        save_pythia_particles(e, pythia_particles, nof_particles, kin);
        e.flag.need_resample_dir = false;
        e.flag.need_resample_phi = false;
      }
      else // SPP from Delta
      {
        // e.flag.res_delta = true;      // mark pion as comming from Delta
        gen_final_particles_hybrid(pion_pdg); // generate final state particles
      }
    } else // more inelastic final state or single kaon production
    {
      e.weight = xsec.get_total(kin.jacobian);
      save_pythia_particles(e, pythia_particles, nof_particles, kin);
      e.flag.need_resample_dir = false;
      e.flag.need_resample_phi = false;
    }
    delete pythia71;
  }
  // set all outgoing particles position to target nucleon position
  for (int j = 0; j < e.out.size(); j++)
    e.out[j].r = e.in[1].r;
}

double hybrid_dsdQ2dW_tab(res_kinematics *kin, int params[4], vect final_pion, double pion_momentum)
{
  double result = 0.;
  double Q2 =-kin->q*kin->q;
  double W  = kin->W;

  // double *hybrid_grid[] = {hybrid_grid_dQ2dW_11,   hybrid_grid_dQ2dW_22,   hybrid_grid_dQ2dW_21,
  //                          hybrid_grid_dQ2dW_2_11, hybrid_grid_dQ2dW_2_22, hybrid_grid_dQ2dW_2_21};
  auto& hybrid_grid = hybrid_grid_mmapio::get2d();
  int hg_idx = hybrid_grid_idx(params[2], params[3]*10 + params[1]);

  double Q2min,Q2max,Q2spc,Q2bin;
  double  Wmin, Wmax, Wspc, Wbin;

  Q2min = Q2min_hybrid; Q2max = Q2max_hybrid;
  Q2spc = Q2spc_hybrid; Q2bin = Q2bin_hybrid;
  Wmin  =  Wmin_hybrid;  Wmax =  Wmax_hybrid;
  Wspc  =  Wspc_hybrid;  Wbin =  Wbin_hybrid;

  if(Q2 < Q2max_2_hybrid) // switch to the denser mesh
  {
    hg_idx += 3;
    Q2min = Q2min_2_hybrid; Q2max = Q2max_2_hybrid;
    Q2spc = Q2spc_2_hybrid; Q2bin = Q2bin_2_hybrid;
  }

  auto l = get_lepton_vec(kin->neutrino.E(), Q2, W, kin->lepton_mass, res_kinematics::avg_nucleon_mass, params);

  // interpolate the nuclear tensor elements
  double w[5] = {0,0,0,0,0}; // 00, 03, 33, 11/22, 12

  // if the variables are within the grid
  if( Q2 >= Q2min_hybrid && Q2 <= Q2max_hybrid && W >= Wmin_hybrid && W <= Wmax_hybrid )
  {
    // bilinear interpolation in axis: x(Q2), y(W), each field is 5 numbers
    int    Q2f = int((Q2-Q2min)/Q2spc); // number of bin in Q2 (floor)
    double Q2d = Q2-Q2min-Q2f*Q2spc;    // distance from the previous point
           Q2d/= Q2spc;                 // normalized
    int     Wf = int((W-Wmin)/Wspc);    // number of bin in W (floor)
    double  Wd = W-Wmin-Wf*Wspc;        // distance from the prefious point
            Wd/= Wspc;                  // normalized
  
    // 9 points surrounding the desired point
    if (Q2f == 0 || (Q2d > 0.5 && Q2f < Q2bin - 2)) {
      Q2f++;
      Q2d--;
    }
    if (Wf == 0 || (Wd > 0.5 && Wf < Wbin - 2)) {
      Wf++;
      Wd--;
    }
    std::array<std::array<double, 3>, 3> vars;
    for (int i = 0; i < 5; i++) {
      for (int k = 0; k < 3; k++)
        for (int l = 0; l < 3; l++) {
          auto index =
              index_calculator_2d(Wf + (k - 1), Q2f + (l - 1), Q2bin) * 5;
          vars[k][l] = hybrid_grid[hg_idx][index + i];
        }
      w[i] = poly_interop_2d(vars, Wd, Q2d);
    }
  }

  // contract the tensors
  result = l[0]*w[0] + 2*l[1]*w[1] + l[2]*w[2] + 0.5*l[3]*w[3] + 2*params[2]*l[4]*w[4];

  return result;
}

double hybrid_dsdQ2dWdcth_tab(res_kinematics *kin, int params[4], vect final_pion, double pion_momentum)
{
  // double result_a = hybrid_dsdQ2dWdcth(kin, params, final_pion, pion_momentum);
  double result = 0.;
  double Q2 =-kin->q*kin->q;
  double W  = kin->W;

  // double *hybrid_grid[] = {hybrid_grid_dQ2dWdcth_11,   hybrid_grid_dQ2dWdcth_22,   hybrid_grid_dQ2dWdcth_21,
  //                          hybrid_grid_dQ2dWdcth_2_11, hybrid_grid_dQ2dWdcth_2_22, hybrid_grid_dQ2dWdcth_2_21};
  auto & hybrid_grid = hybrid_grid_mmapio::get();
  int hg_idx = hybrid_grid_idx(params[2], params[3]*10 + params[1]);

  double  Q2min,  Q2max,  Q2spc,  Q2bin;
  double   Wmin,   Wmax,   Wspc,   Wbin;
  double cthmin, cthmax, cthspc, cthbin;

  Q2min  =  Q2min_hybrid; Q2max  =  Q2max_hybrid;
  Q2spc  =  Q2spc_hybrid; Q2bin  =  Q2bin_hybrid;
  Wmin   =   Wmin_hybrid; Wmax   =   Wmax_hybrid;
  Wspc   =   Wspc_hybrid; Wbin   =   Wbin_hybrid;
  cthmin = cthmin_hybrid; cthmax = cthmax_hybrid;
  cthspc = cthspc_hybrid; cthbin = cthbin_hybrid;

  if(Q2 < Q2max_2_hybrid) // switch to the denser mesh
  {
    hg_idx += 3;
    Q2min = Q2min_2_hybrid; Q2max = Q2max_2_hybrid;
    Q2spc = Q2spc_2_hybrid; Q2bin = Q2bin_2_hybrid;
  }

  // find proper angle costh_pi^ast
  vect k = kin->neutrino;  //
  vect kp= kin->lepton;    // They are in target rest frame!
  vect q = kin->q;         //
  k.boost (-kin->hadron_speed);
  kp.boost(-kin->hadron_speed);
  q.boost (-kin->hadron_speed);
  vec Zast = q;
  vec Yast = vecprod(k,kp);
  vec Xast = vecprod(Yast,q);
  Zast.normalize(); Yast.normalize(); Xast.normalize();
  double pion_cos_theta = Zast * vec(final_pion) / final_pion.length();

  auto l = get_lepton_vec(kin->neutrino.E(), Q2, W, kin->lepton_mass, res_kinematics::avg_nucleon_mass, params);

  // interpolate the nuclear tensor elements
  double w[5] = {0,0,0,0,0}; // 00, 03, 33, 11/22, 12

  // if the variables are within the grid
  if( Q2 >= Q2min_hybrid && Q2 <= Q2max_hybrid && W >= Wmin_hybrid && W <= Wmax_hybrid &&
      pion_cos_theta >= cthmin_hybrid && pion_cos_theta <= cthmax_hybrid )
  {
    // trilinear interpolation in axis: x(Q2), y(W), z(costh), each field is 5 numbers
    int     Q2f = int((Q2-Q2min)/Q2spc);               // number of bin in Q2 (floor)
    double  Q2d = Q2-Q2min-Q2f*Q2spc;                  // distance from the previous point
            Q2d/= Q2spc;                               // normalized
    int      Wf = int((W-Wmin)/Wspc);                  // number of bin in W (floor)
    double   Wd = W-Wmin-Wf*Wspc;                      // distance from the previous point
             Wd/= Wspc;                                // normalized
    int    cthf = int((pion_cos_theta-cthmin)/cthspc); // number of bin in cth (floor)
    double cthd = pion_cos_theta-cthmin-cthf*cthspc;   // distance from the previous point
           cthd/= cthspc;                              // normalized
    // normalize distance to [-0.5, 0.5], modify the index correspondingly
    // Q2d, Wd, cthd as 0-point, to be used for interpolation
    if (Q2f==0 || (Q2d > 0.5 && Q2f < Q2bin-2)){
      Q2f++;
      Q2d--;
    }
    if (Wf == 0 || (Wd > 0.5 && Wf < Wbin - 2)) {
      Wf++;
      Wd--;
    }
    if (cthf == 0 || (cthd > 0.5 && cthf < cthbin - 2)) {
      cthf++;
      cthd--;
    }
    std::array<std::array<std::array<double, 3>, 3>, 3> vars;
    for (int i = 0; i < 5; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++) {
            auto index = index_calculator(cthf + (j - 1), Wf + (k - 1),
                                          Q2f + (l - 1), Q2bin) *
                         5;
            vars[j][k][l] = hybrid_grid[hg_idx][index + i];
          }
        }
      }
      w[i] = poly_interop_3d(vars, cthd, Wd, Q2d);
    }
  }
  else {
    __builtin_trap();
  }

  // contract the tensors
  result = l[0]*w[0] + 2*l[1]*w[1] + l[2]*w[2] + 0.5*l[3]*w[3] + 2*params[2]*l[4]*w[4]; // *2Pi

  result *= pion_momentum / pow(2*Pi,3);                                                // /2Pi
  result *= 2; // Phase space!
  // if (abs(result_a - result) / (result_a + result) > 0.05){
  //   std::cerr << "result_a = " << result_a << " result = " << result << std::endl;
  //   for ( int i = 0; i < 5; i++ ) std::cerr << "w[" << i << "] = " << w[i] << std::endl;
  //   for ( int i = 0; i < 5; i++ ) std::cerr << "w_global[" << i << "] = " << w_global[i] << std::endl;
  //   for ( int i = 0; i < 5; i++ ) std::cerr << "dw[" << i << "] = " << (w[i] - w_global[i])/(w[i] + w_global[i]) << std::endl;
  //   for ( int i = 0; i < 5; i++ ) std::cerr << "l[" << i << "] = " << l[i] << std::endl;
  //   for ( int i = 0; i < 5; i++ ) std::cerr << "l_global[" << i << "] = " << l_global[i] << std::endl;
  //   for ( int i = 0; i < 5; i++ ) std::cerr << "dl[" << i << "] = " << (l[i] - l_global[i])/(l[i] + l_global[i]) << std::endl;
  //   raise(SIGTRAP);
  // }
  return result;
}

double hybrid_dsdQ2dWdcth(res_kinematics* kin, int params[4], vect final_pion, double pion_momentum)
{
  // placeholders
  double costh[1];
  double ABCDE[1][5] = {{0,0,0,0,0}};
  double result = 0.;

  // get Q^2, W
  double Q2 =-kin->q*kin->q;
  double W  = kin->W;

  // find proper angles costh_pi^ast and phi_pi^ast
  vect k = kin->neutrino;  //
  vect kp= kin->lepton;    // They are in target rest frame!
  vect q = kin->q;         //
  k.boost (-kin->hadron_speed);
  kp.boost(-kin->hadron_speed);
  q.boost (-kin->hadron_speed);
  vec Zast = q;
  vec Yast = vecprod(k,kp);
  vec Xast = vecprod(Yast,q);
  Zast.normalize(); Yast.normalize(); Xast.normalize();
  double pion_cos_theta = Zast * vec(final_pion) / final_pion.length();
  double pion_phi = atan2(Yast*vec(final_pion),Xast*vec(final_pion));

  // fill costh
  costh[0] = pion_cos_theta;

  // get ABCDE
  hybrid_ABCDE(kin->neutrino.E(), Q2, W, kin->lepton_mass, res_kinematics::avg_nucleon_mass, costh, 1, params, ABCDE, ABCDE);

  result  = ABCDE[0][0];                 // *2Pi
  result *= pion_momentum / pow(2*Pi,3); // /2Pi
  result *= 2; // Phase space!

  return result;
}

double hybrid_dsdQ2dWdOm(res_kinematics* kin, int params[4], vect final_pion, double pion_momentum)
{
  // placeholders
  double costh[1];
  double ABCDE[1][5] = {{0,0,0,0,0}};
  double result = 0.;

  // get Q^2, W
  double Q2 =-kin->q*kin->q;
  double W  = kin->W;

  // find proper angles costh_pi^ast and phi_pi^ast
  vect k = kin->neutrino;  //
  vect kp= kin->lepton;    // They are in target rest frame!
  vect q = kin->q;         //
  k.boost (-kin->hadron_speed);
  kp.boost(-kin->hadron_speed);
  q.boost (-kin->hadron_speed);
  vec Zast = q;
  vec Yast = vecprod(k,kp);
  vec Xast = vecprod(Yast,q);
  Zast.normalize(); Yast.normalize(); Xast.normalize();
  double pion_cos_theta = Zast * vec(final_pion) / final_pion.length();
  double pion_phi = atan2(Yast*vec(final_pion),Xast*vec(final_pion));

  // fill costh
  costh[0] = pion_cos_theta;

  // get ABCDE
  hybrid_ABCDE(kin->neutrino.E(), Q2, W, kin->lepton_mass, res_kinematics::avg_nucleon_mass, costh, 1, params, ABCDE, ABCDE);

  result  = ABCDE[0][0] + ABCDE[0][1]*cos(pion_phi) + ABCDE[0][2]*cos(2*pion_phi)
                        + ABCDE[0][3]*sin(pion_phi) + ABCDE[0][4]*sin(2*pion_phi);
  result *= pion_momentum / pow(2*Pi,4);
  result *= 4 * Pi; // Phase space!

  return result;
}

double hybrid_sample_costh(double Enu, double Q2, double W, double m, double M, int params[4], int costh_pts)
{
  // Choosing cos_th^* from dsdQ2dWdcosth
  double costh_rnd;

  // Specify points for the polynomial interpolation
  //const int costh_mem = 3;             // maximum number of points: 3, 5, 7, 9, ...
  //int costh_pts = 3;                   // actual number of points:  3, 5, 7, 9, ...
  //double costh[costh_mem];             // Points of interpolation

  // Allocate memory dynamically
  double *costh = new double[costh_pts];

  for( int i = 0; i < costh_pts; i++ )           // Fill costh with evenly spaced points
    costh[i] = -0.75 + i * 1.5/(costh_pts-1);    // from -0.75 to 0.75
    //costh[i] = -cos( Pi / (costh_pts-1) * i ); // from -1 to 1

  // Get the A function, \propto dsdQ2dWdcosth
  //double ABCDE[costh_mem][5]; double ds[costh_mem];
  double (*ABCDE)[5] = new double[costh_pts][5];
  double *ds = new double[costh_pts];
  hybrid_ABCDE(Enu, Q2, W, m, M, costh, costh_pts, params, ABCDE, ABCDE);
  for( int i = 0; i < costh_pts; i++ )
    ds[i] = ABCDE[i][0];

  // Fit a polynomial to given number of points in ds(costh)
  //double poly_coeffs[costh_mem];
  double *poly_coeffs = new double[costh_pts];
  hybrid_poly_fit(costh_pts, costh, ds, poly_coeffs);     // use Root methods
  //hybrid_poly_fit_2(costh_pts, costh, ds, poly_coeffs); // analytical

  // Normalize the coefficients to obtain a probability density
  double norm = hybrid_poly_dist(costh_pts, poly_coeffs, cthmin_hybrid, cthmax_hybrid);
  for( int i = 0; i < costh_pts; i++ )
    poly_coeffs[i] /= norm;

  // Choose a value for costh_rnd
  if( costh_pts == 3 ) // analytical
    costh_rnd = hybrid_poly_rnd_2(costh_pts, poly_coeffs, cthmin_hybrid, cthmax_hybrid);
  if( costh_pts > 3 || fabs(costh_rnd) > 1 ) // bisection
    costh_rnd = hybrid_poly_rnd(costh_pts, poly_coeffs, cthmin_hybrid, cthmax_hybrid, 0.00001);

  // Free dynamically allocated memory
  delete [] costh;
  delete [] ABCDE;
  delete [] ds;
  delete [] poly_coeffs;

  return costh_rnd;
}

double hybrid_sample_costh_2(double Enu, double Q2, double W, double M, int params[4], vect kl_inc_lab, vect kl_lab)
{
  // Choosing cos_th^* from dsdQ2dWdcosth
  double costh_rnd;

  // Specify all details needed for arrays
  // double *hybrid_grid[] = {hybrid_grid_dQ2dWdcth_11,   hybrid_grid_dQ2dWdcth_22,   hybrid_grid_dQ2dWdcth_21,
  //                          hybrid_grid_dQ2dWdcth_2_11, hybrid_grid_dQ2dWdcth_2_22, hybrid_grid_dQ2dWdcth_2_21};
  auto & hybrid_grid = hybrid_grid_mmapio::get();
  int hg_idx = hybrid_grid_idx(params[2], params[3]*10 + params[1]);

  double  Q2min,  Q2max,  Q2spc,  Q2bin;
  double   Wmin,   Wmax,   Wspc,   Wbin;
  double cthmin, cthmax, cthspc, cthbin;

  Q2min  =  Q2min_hybrid; Q2max  =  Q2max_hybrid;
  Q2spc  =  Q2spc_hybrid; Q2bin  =  Q2bin_hybrid;
  Wmin   =   Wmin_hybrid; Wmax   =   Wmax_hybrid;
  Wspc   =   Wspc_hybrid; Wbin   =   Wbin_hybrid;
  cthmin = cthmin_hybrid; cthmax = cthmax_hybrid;
  cthspc = cthspc_hybrid; cthbin = cthbin_hybrid;

  if(Q2 < Q2max_2_hybrid) // switch to the denser mesh
  {
    hg_idx += 3;
    Q2min = Q2min_2_hybrid; Q2max = Q2max_2_hybrid;
    Q2spc = Q2spc_2_hybrid; Q2bin = Q2bin_2_hybrid;
  }

  // we build leptonic tensor in CMS with q along z, see JES paper App. A
  vect q_lab  = kl_inc_lab - kl_lab;       // They are not in LAB, but in target rest frame!

  double v = q_lab.length() / (q_lab[0] + M);
  double g = 1 / sqrt(1 - v*v);
  double c = (kl_inc_lab.length() - kl_lab[3]) / q_lab.length();
  double s = sqrt(pow(kl_lab.length(),2) - pow(kl_lab[3],2)) / q_lab.length();

  vect kl_inc (g*(kl_inc_lab[0] - v*kl_inc_lab.length()*c), kl_inc_lab.length()*s,
               0, g*(-v*kl_inc_lab[0] + kl_inc_lab.length()*c));
  vect kl     (g*(kl_lab[0] - v*(kl_inc_lab.length()*c - q_lab.length())), kl_inc_lab.length()*s,
               0, g*(-v*kl_lab[0] + kl_inc_lab.length()*c - q_lab.length()));

  double kl_inc_dot_kl = kl_inc * kl;

  // A function placeholder
  double *A = new double[cthbin_hybrid];

  // calculate the leptonic tensor elements
  double l[5] = {0,0,0,0,0}; // 00, 03, 33, 11+22, 12
  l[0] = (2*kl_inc[0]*kl[0] - kl_inc_dot_kl);
  l[1] = ( -kl_inc[0]*kl[3] - kl[0]*kl_inc[3]);
  l[2] = (2*kl_inc[3]*kl[3] + kl_inc_dot_kl);
  l[3] = (2*kl_inc[1]*kl[1] + kl_inc_dot_kl);
  l[3]+= (2*kl_inc[2]*kl[2] + kl_inc_dot_kl);
  l[4] = (  kl_inc[0]*kl[3] - kl[0]*kl_inc[3]);

  // interpolate the nuclear tensor elements
  double w[5] = {0,0,0,0,0}; // 00, 03, 33, 11/22, 12

  // if the variables are within the grid
  if( Q2 >= Q2min_hybrid && Q2 <= Q2max_hybrid && W >= Wmin_hybrid && W <= Wmax_hybrid )
  {
    // trilinear interpolation in axis: x(Q2), y(W), z(costh), each field is 5 numbers
    int     Q2f = int((Q2-Q2min)/Q2spc);               // number of bin in Q2 (floor)
    double  Q2d = Q2-Q2min-Q2f*Q2spc;                  // distance from the previous point
            Q2d/= Q2spc;                               // normalized
    int      Wf = int((W-Wmin)/Wspc);                  // number of bin in W (floor)
    double   Wd = W-Wmin-Wf*Wspc;                      // distance from the previous point
             Wd/= Wspc;                                // normalized

    for( int cthf = 0; cthf < cthbin_hybrid; cthf++ ) // for all cth points
    {
      // 4 points surrounding the desired point
      int p00 = (cthf*Wbin*Q2bin+Wf*Q2bin+Q2f)*5; // bottom left
      int p10 = p00+5;                            // bottom right
      int p01 = p00+Q2bin*5;                      // top left
      int p11 = p01+5;                            // top right

      // interpolate
      for( int i = 0; i < 5; i++ )
        w[i] = bilinear_interp(hybrid_grid[hg_idx][p00+i], hybrid_grid[hg_idx][p10+i],
                               hybrid_grid[hg_idx][p01+i], hybrid_grid[hg_idx][p11+i], Q2d, Wd);
      // contract the tensors
      A[cthf] = l[0]*w[0] + 2*l[1]*w[1] + l[2]*w[2] + 0.5*l[3]*w[3] + 2*params[2]*l[4]*w[4];
    }
  }

  // find maximum
  double maximum = -1;
  for( int i = 0; i < cthbin_hybrid; i++ )
  {
    if( A[i] > maximum )
      maximum = A[i];
  }

  // accept-or-reject
  double A_value;
  do
  {
      costh_rnd = frandom()*(cthmax-cthmin)+cthmin;
    int    cthf = int((costh_rnd-cthmin)/cthspc); // number of bin in cth (floor)
    double cthd = costh_rnd-cthmin-cthf*cthspc;   // distance from the previous point
           cthd/= cthspc;                         // normalized
        A_value = linear_interp(A[cthf], A[cthf+1], cthd);
  }
  while( A_value/maximum < frandom() );

  // Free dynamically allocated memory
  delete [] A;

  return costh_rnd;
}

double hybrid_sample_phi(double Enu, double Q2, double W, double m, double M, int params[4], double costh_rnd)
{
  // Choosing phi^* from dsdQ2dWdOm
  double phi_rnd;

  // Specify points for the polynomial interpolation
  const int phi_pts = 3;             // 3, 5, 7, 9, ...
  double phi[phi_pts];               // Points of interpolation
  for( int i = 0; i < phi_pts; i++ ) // Fill phi with evenly spaced points
    phi[i] = 2*Pi / (phi_pts-1) * i - Pi;

  // Get the ABCDE combination, \propto dsdQ2dWdOm
  double ABCDE[1][5]; double ds[phi_pts];
  hybrid_ABCDE(Enu, Q2, W, m, M, &costh_rnd, 1, params, ABCDE, ABCDE);
  for( int i = 0; i < phi_pts; i++ )
    ds[i] = ABCDE[0][0] + ABCDE[0][1]*cos(phi[i]) + ABCDE[0][2]*cos(2*phi[i])
                        + ABCDE[0][3]*sin(phi[i]) + ABCDE[0][4]*sin(2*phi[i]);

  // Fit a polynomial to given number of points in ds(phi)
  double poly_coeffs[phi_pts];
  hybrid_poly_fit(phi_pts, phi, ds, poly_coeffs);

  // Normalize the coefficients to obtain a probability density
  double norm = hybrid_poly_dist(phi_pts, poly_coeffs, phi[0], phi[phi_pts-1]);
  for( int i = 0; i < phi_pts; i++ )
    poly_coeffs[i] /= norm;

  // Choose a value for phi_rnd
  phi_rnd = hybrid_poly_rnd(phi_pts, poly_coeffs, phi[0], phi[phi_pts-1], 0.001);

  return phi_rnd;
}

double hybrid_sample_phi_2(double Enu, double Q2, double W, double m, double M, int params[4], double costh_rnd)
{
  // Choosing phi^* from dsdQ2dWdOm
  double phi_rnd;

  // Get the ABCDE combination, \propto dsdQ2dWdOm
  double ABCDE[1][5];
  hybrid_ABCDE(Enu, Q2, W, m, M, &costh_rnd, 1, params, ABCDE, ABCDE);

  // Integral is
  // A*x + B*sinx + 0.5*C*sin2x - D*cosx - 0.5*E*cos2x

  // Normalize the coefficients to obtain a probability density
  double norm = ABCDE[0][0]*2*Pi;
  for( int i = 0; i < 5; i++ )
    ABCDE[0][i] /= norm;

  // Choose a value for phi_rnd
  phi_rnd = hybrid_dcmp_rnd(ABCDE, -Pi, Pi, 0.001);

  return phi_rnd;
}

double hybrid_sample_phi_3(double Enu, double Q2, double W, double m, double M, int params[4], double costh_rnd)
{
  // Choosing phi^* from dsdQ2dWdOm
  double phi_rnd;

  // Get the ABCDE combination, \propto dsdQ2dWdOm
  double ABCDE[1][5];
  hybrid_ABCDE(Enu, Q2, W, m, M, &costh_rnd, 1, params, ABCDE, ABCDE);

  // Integral is
  // A*x + B*sinx + 0.5*C*sin2x - D*cosx - 0.5*E*cos2x

  // Normalize the coefficients to obtain a probability density
  double norm = ABCDE[0][0]*2*Pi;
  for( int i = 0; i < 5; i++ )
    ABCDE[0][i] /= norm;

  // Choose a value for phi_rnd
  phi_rnd = hybrid_dcmp_rnd_2(ABCDE, -Pi, Pi, 0.001);

  return phi_rnd;
}

void hybrid_poly_fit(const int N, double* xpts, double* ypts, double* coeffs)
{
  // We perform a polynomial interpolation
  // Inspired by https://en.wikipedia.org/wiki/Polynomial_interpolation

  // Create matrix of x powers
  TMatrixD A(N,N);
  for( int i = 0; i < N; i++ ) // rows
  {
    for( int j = 0; j < N; j++ ) // columns backwards
    {
      A[i][j] = pow(xpts[i],(N-1-j));
    }
  }

  // Perform an LU decomposition
  TDecompLU LU(A);

  // Create a vector of y points
  TVectorD b(N);
  for( int i = 0; i < N; i++ )
    b[i] = ypts[i];

  // Solve the equation
  LU.Solve(b);
  for( int i = 0; i < N; i ++ )
    coeffs[i] = b[i];
}

void hybrid_poly_fit_2(const int N, double* xpts, double* ypts, double* coeffs)
{
  assert( N == 3 );

  // Calculate the constants in A*x^2 + B*x + C = 0 going through 3 points
  // whose x,y values are in the input arrays

  double denom = (xpts[0] - xpts[1]) * (xpts[0] - xpts[2]) * (xpts[1] - xpts[2]);
  coeffs[0]    = (xpts[2] * (ypts[1] - ypts[0])
                + xpts[1] * (ypts[0] - ypts[2])
                + xpts[0] * (ypts[2] - ypts[1])) / denom;
  coeffs[1]    = (xpts[2] * xpts[2] * (ypts[0] - ypts[1])
                + xpts[1] * xpts[1] * (ypts[2] - ypts[0])
                + xpts[0] * xpts[0] * (ypts[1] - ypts[2])) / denom;
  coeffs[2]    = (xpts[1] * xpts[2] * (xpts[1] - xpts[2]) * ypts[0]
                + xpts[2] * xpts[0] * (xpts[2] - xpts[0]) * ypts[1]
                + xpts[0] * xpts[1] * (xpts[0] - xpts[1]) * ypts[2]) / denom;
}

double hybrid_poly_dist(const int N, double* coeffs, double x_min, double x)
{
  double result = 0.;
  for( int i = 0; i < N; i++ )
  {
    result += coeffs[i]/(N-i) * (pow(x,N-i)-pow(x_min,N-i));
  }
  return result;
}

double hybrid_poly_rnd(const int N, double* coeffs, double x_min, double x_max, double epsilon)
{
  double a = x_min, b = x_max;
  double x, fx;
  double y = frandom();

  do
  {
    x = (b-a)/2 + a;
    fx = hybrid_poly_dist(N, coeffs, x_min, x);
    if( fx > y )
      b = x;
    else
      a = x;
  }
  while( fabs(fx - y) > epsilon );

  return x;
}

double hybrid_poly_rnd_2(const int N, double* coeffs, double x_min, double x_max)
{
  assert( N == 3 );

  double a = coeffs[0] / 3;
  double b = coeffs[1] / 2;
  double c = coeffs[2];
  double d = a - b + c;
  d -= frandom();

  double sol;

  // Solving a cubic equation using the discriminant approach
  // We only need real solutions, no complex parts are calculated.
  // Moreover, the real solution should be between -1 and 1 (this is not garanteed!)
  b /= a;
  c /= a;
  d /= a;

  double disc, q, r, dum1, s, t, term1, r13;
  q  = (3.0*c - (b*b))/9.0;
  r  = -(27.0*d) + b*(9.0*c - 2.0*(b*b));
  r /= 54.0;
  disc  = q*q*q + r*r;
  term1 = (b/3.0);

  double x1_real, x2_real, x3_real;
  if (disc > 0)   // One root real, two are complex
  {
    s = r + sqrt(disc);
    s = s<0 ? -cbrt(-s) : cbrt(s);
    t = r - sqrt(disc);
    t = t<0 ? -cbrt(-t) : cbrt(t);
    sol = -term1 + s + t;
  }
  // The remaining options are all real
  else if (disc == 0)  // All roots real, at least two are equal.
  {
    r13 = r<0 ? -cbrt(-r) : cbrt(r);
    x1_real = -term1 + 2.0*r13;
    if (abs(x1_real) <= 1)
    {
      sol = x1_real;
    }else{
      sol = -(r13 + term1);
    }
  }
  // Only option left is that all roots are real and unequal (to get here, q < 0)
  else
  {
    q = -q;
    dum1 = q*q*q;
    dum1 = acos(r/sqrt(dum1));
    r13 = 2.0*sqrt(q);
    x1_real = -term1 + r13*cos(dum1/3.0);
    x2_real = -term1 + r13*cos((dum1 + 2.0*Pi)/3.0);
    if (abs(x1_real) <= 1.){sol = x1_real;}
    else if( abs(x2_real) <= 1){sol = x2_real;}
    else{sol = -term1 + r13*cos((dum1 + 4.0*Pi)/3.0);}
  }

  return sol;
}

double hybrid_dcmp_rnd(double (*ABCDE)[5], double x_min, double x_max, double epsilon)
{
  double a = x_min, b = x_max;
  double x, fx;
  double y = frandom();

  do
  {
    x = (b-a)/2 + a;
    fx = ABCDE[0][0]*(x-x_min) + ABCDE[0][1]*(sin(x)-sin(x_min)) + ABCDE[0][2]*(sin(2*x)-sin(2*x_min))/2
                               - ABCDE[0][3]*(cos(x)-cos(x_min)) - ABCDE[0][4]*(cos(2*x)-cos(2*x_min))/2;
    if( fx > y )
      b = x;
    else
      a = x;
  }
  while( fabs(fx - y) > epsilon );

  return x;
}

double hybrid_dcmp_rnd_2(double (*ABCDE)[5], double x_min, double x_max, double epsilon)
{
  double a = x_min, b = x_max;
  double x, fx, fxp;
  double y = frandom();

  x = (x_max-x_min)/2 + x_min;

  int stp = 0;
  do
  {
    fx = ABCDE[0][0]*(x-x_min) + ABCDE[0][1]*(sin(x)-sin(x_min)) + ABCDE[0][2]*(sin(2*x)-sin(2*x_min))/2
                               - ABCDE[0][3]*(cos(x)-cos(x_min)) - ABCDE[0][4]*(cos(2*x)-cos(2*x_min))/2;
    fxp = ABCDE[0][0] + ABCDE[0][1]*cos(x) + ABCDE[0][2]*cos(2*x)
                      + ABCDE[0][3]*sin(x) + ABCDE[0][4]*sin(2*x);
    double new_x = x - (fx - y)/fxp;
    // if the function is very flat and the next step go out of bounds
    // use one bisective step to get out
    if( fabs(fxp) < epsilon || new_x < a || new_x > b || stp > 9 )
    {
      if( fx > y )
        b = x;
      else
        a = x;
      x = (b-a)/2 + a;
      stp = 0;
    }
    else
    {
      x = new_x;
      stp++;
    }
    //assert(stp<20);
  }
  while( fabs(fx - y) > epsilon );

  return x;
}

vec hybrid_dir_from_adler(double costh, double phi, vect k, vect kp)
{
  vec Zast = k - kp;
  vec Yast = vecprod(k,kp);
  vec Xast = vecprod(Yast,Zast);
  Zast.normalize(); Yast.normalize(); Xast.normalize();

  vec kierunek = costh * Zast + sqrt(1 - costh*costh) * (cos(phi) * Xast + sin(phi) * Yast);

  return kierunek;
}

void resevent_dir_hybrid(event& e, nucleus& t, int method)
{
  // get all 4-vectors from the event and boost them to the N-rest frame
  particle target = e.in[1];  
  target.t       -= get_binding_energy (e.par, target, t);
  vect neutrino   = e.in[0];  neutrino.boost(-target.v());
  vect lepton     = e.out[0]; lepton.boost  (-target.v());
  vect pion       = e.out[1]; pion.boost    (-target.v());
  vect nucleon    = e.out[2]; nucleon.boost (-target.v());

  // calculate hadron_speed and boost to CMS
  vect q = neutrino - lepton;
  double Mef = min(sqrt(target.p4() * target.p4()), res_kinematics::avg_nucleon_mass);
  double q3  = sqrt(pow(Mef + q.t,2) - e.W()*e.W());
  vec hadron_speed = q / sqrt(e.W()*e.W() + q3 * q3);
  pion.boost     (-hadron_speed);
  nucleon.boost  (-hadron_speed);

  // generate params for the Ghent code
  int params[4];
  params[0] = 1;                                   // only CC for now
  params[2] = (1 - 2.0 * (e.out[0].pdg > 0));
  if( e.in[0].pdg > 0 ) // neutrino
  {
    if( e.in[1].pdg == PDG::pdg_proton )
    {
      params[3] = 1; params[1] = 1;
    }
    else
    {
      params[3] = 2;
      if( e.out[1].pdg == PDG::pdg_pi )
        params[1] = 1;
      else
        params[1] = 2;
    }
  }
  else                   // antineutrino
  {
    if( e.in[1].pdg == PDG::pdg_neutron )
    {
      params[3] = 2; params[1] = 2;
    }
    else
    {
      params[3] = 1;
      if( e.out[1].pdg == PDG::pdg_pi )
        params[1] = 1;
      else
        params[1] = 2;
    }
  }

  // Choose cos_theta^* for dsdQ2dW, in Adler frame. E_nu in N-rest!
  double costh_rnd;
  if(e.W() > Wmax_hybrid || -e.q2() > Q2max_hybrid) // outside of limits, dsdQ2dWdcth was used
  {
    vec Zast  = neutrino - lepton;
    vec Yast  = vecprod(neutrino,lepton);
    vec Xast  = vecprod(Yast,Zast);
    Zast.normalize(); Yast.normalize(); Xast.normalize();
    costh_rnd = Zast * vec(pion) / pion.length();
  }
  else
  {
    if( method > 0 )
      costh_rnd = hybrid_sample_costh(neutrino.t, -e.q2(), e.W(), e.out[0].m(), Mef, params, method);
    else
      costh_rnd = hybrid_sample_costh_2(neutrino.t, -e.q2(), e.W(), Mef, params, neutrino, lepton);
  }

  // Choose phi^* for dsdQ2dWdcosth, in Adler frame. E_nu in N-rest!
  //double phi_rnd = hybrid_sample_phi(neutrino.t, -e.q2(), e.W(), e.out[0].m(), Mef, params, costh_rnd);
  //double phi_rnd = hybrid_sample_phi_2(neutrino.t, -e.q2(), e.W(), e.out[0].m(), Mef, params, costh_rnd);
  double phi_rnd = hybrid_sample_phi_3(neutrino.t, -e.q2(), e.W(), e.out[0].m(), Mef, params, costh_rnd);

  // Modify final hadron directions as specified in the Adler frame
  neutrino.boost (-hadron_speed); // Neutrino and lepton have to be boosted to CMS
  lepton.boost   (-hadron_speed); // This is needed to properly specify the Adler frame
  vec dir_rnd = hybrid_dir_from_adler(costh_rnd, phi_rnd, neutrino, lepton);

  // Recalculate the hadronic kinematics, dir_rnd is the new direction of pion
  double momentum = nucleon.length();
  nucleon = vect(nucleon.t, -momentum * dir_rnd.x, -momentum * dir_rnd.y, -momentum * dir_rnd.z);
  pion = vect(pion.t, momentum * dir_rnd.x, momentum * dir_rnd.y, momentum * dir_rnd.z);

  // Boost hadrons back
  nucleon.boost(hadron_speed);
  nucleon.boost(target.v());
  pion.boost   (hadron_speed);
  pion.boost   (target.v());

  // Correct the particles in the out vector
  e.out[1].p4() = pion;
  e.out[2].p4() = nucleon;
}

void resevent_phi_hybrid(event& e, nucleus& t)
{
  // get all 4-vectors from the event and boost them to the N-rest frame
  particle target = e.in[1];  
  target.t       -= get_binding_energy (e.par, target, t);
  vect neutrino   = e.in[0];  neutrino.boost(-target.v());
  vect lepton     = e.out[0]; lepton.boost  (-target.v());
  vect pion       = e.out[1]; pion.boost    (-target.v());
  vect nucleon    = e.out[2]; nucleon.boost (-target.v());

  // calculate hadron_speed and boost to CMS
  vect q = neutrino - lepton;
  double Mef = min(sqrt(target.p4() * target.p4()), res_kinematics::avg_nucleon_mass);
  double q3  = sqrt(pow(Mef + q.t,2) - e.W()*e.W());
  vec hadron_speed = q / sqrt(e.W()*e.W() + q3 * q3);
  pion.boost     (-hadron_speed);
  nucleon.boost  (-hadron_speed);

  // generate params for the Ghent code
  int params[4];
  params[0] = 1;                                   // only CC for now
  params[2] = (1 - 2.0 * (e.out[0].pdg > 0));
  if( e.in[0].pdg > 0 ) // neutrino
  {
    if( e.in[1].pdg == PDG::pdg_proton )
    {
      params[3] = 1; params[1] = 1;
    }
    else
    {
      params[3] = 2;
      if( e.out[1].pdg == PDG::pdg_pi )
        params[1] = 1;
      else
        params[1] = 2;
    }
  }
  else                   // antineutrino
  {
    if( e.in[1].pdg == PDG::pdg_neutron )
    {
      params[3] = 2; params[1] = 2;
    }
    else
    {
      params[3] = 1;
      if( e.out[1].pdg == PDG::pdg_pi )
        params[1] = 1;
      else
        params[1] = 2;
    }
  }

  // Choose cos_theta^* for dsdQ2dW, in Adler frame. E_nu in N-rest!
  //double costh_rnd = hybrid_sample_costh(neutrino.t, -e.q2(), e.W(), e.out[0].m(), params);
  vec Zast = neutrino - lepton;
  vec Yast = vecprod(neutrino,lepton);
  vec Xast = vecprod(Yast,Zast);
  Zast.normalize(); Yast.normalize(); Xast.normalize();
  double costh = Zast * vec(pion) / pion.length();

  // Choose phi^* for dsdQ2dWdcosth, in Adler frame. E_nu in N-rest!
  //double phi_rnd = hybrid_sample_phi(neutrino.t, -e.q2(), e.W(), e.out[0].m(), Mef, params, costh_rnd);
  //double phi_rnd = hybrid_sample_phi_2(neutrino.t, -e.q2(), e.W(), e.out[0].m(), Mef, params, costh_rnd);
  double phi_rnd = hybrid_sample_phi_3(neutrino.t, -e.q2(), e.W(), e.out[0].m(), Mef, params, costh);

  // Modify final hadron directions as specified in the Adler frame
  neutrino.boost (-hadron_speed); // Neutrino and lepton have to be boosted to CMS
  lepton.boost   (-hadron_speed); // This is needed to properly specify the Adler frame
  vec dir_rnd = hybrid_dir_from_adler(costh, phi_rnd, neutrino, lepton);

  // Recalculate the hadronic kinematics, dir_rnd is the new direction of pion
  double momentum = nucleon.length();
  nucleon = vect(nucleon.t, -momentum * dir_rnd.x, -momentum * dir_rnd.y, -momentum * dir_rnd.z);
  pion = vect(pion.t, momentum * dir_rnd.x, momentum * dir_rnd.y, momentum * dir_rnd.z);

  // Boost hadrons back
  nucleon.boost(hadron_speed);
  nucleon.boost(target.v());
  pion.boost   (hadron_speed);
  pion.boost   (target.v());

  // Correct the particles in the out vector
  e.out[1].p4() = pion;
  e.out[2].p4() = nucleon;
}

double linear_interp(double f0, double f1, double xd)
{
  return f0 * (1-xd) + f1 * xd;
}

double bilinear_interp(double f00,  double f10,  double f01,  double f11,  double xd, double yd)
{
  return linear_interp(linear_interp(f00, f10, xd), linear_interp(f01, f11, xd), yd);
}

double trilinear_interp(double f000, double f010, double f100, double f110,
                        double f001, double f011, double f101, double f111, double xd, double yd, double zd)
{
  return linear_interp(bilinear_interp(f000, f010, f100, f110, xd, yd),
                       bilinear_interp(f001, f011, f101, f111, xd ,yd), zd);
}

int hybrid_grid_idx(int helicity, int channel)
{
  int hg_idx = -1;

  if(helicity < 0)
  {
    switch (channel)
    {
      case 11: hg_idx = 0; break;
      case 22: hg_idx = 1; break;
      case 21: hg_idx = 2; break;
      default:
        cerr << "[WARNING]: no hybrid inclusive grid!\n";
    }
  }
  else
  {
    switch (channel)
    {
      case 11: hg_idx = 2; break; // anu_11 has the tables of nu_21
      case 12: hg_idx = 1; break; // anu_12 has the tables of nu_22
      case 22: hg_idx = 0; break; // anu_21 has the tables of nu_11
      default:
        cerr << "[WARNING]: no hybrid inclusive grid!\n";
    }
  }

  return hg_idx;
}
