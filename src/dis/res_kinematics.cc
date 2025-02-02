#include "res_kinematics.h"
#include "jednostki.h"
#include "dis/LeptonMass.h"
#include "nucleus.h"

const double res_kinematics::Wmin = 1080;  // TODO: it is not exactly pion mass + nucleon mass
const double res_kinematics::avg_nucleon_mass = (PDG::mass_proton + PDG::mass_neutron) / 2.0;
const double res_kinematics::pythia_threshold = 1210;

inline double pow2(double x) { return x * x; }

// set all necessary variables so is_above_threshold may be called
res_kinematics::res_kinematics(const event &e, nucleus &t) : neutrino(e.in[0]), target(e.in[1]) {
  // final lepton mass = 0 for NC or corresponding lepton mass (nu PDG - 1)
  lepton_mass = e.flag.cc * PDG::mass(abs(neutrino.pdg) - 1);
  lepton_mass2 = lepton_mass * lepton_mass;

  // subtract binding energy from nucleon energy inside nucleus
  target.t -= get_binding_energy(e.par, target, t);

  // boost to the bound nucleon rest frame
  neutrino.boost(-target.v());

  // effective nucleon mass depends on binding energy (cannot be larger than avg mass of nucleon)
  effective_mass = min(sqrt(target.p4() * target.p4()), avg_nucleon_mass);
  effective_mass2 = effective_mass * effective_mass;
}

bool res_kinematics::generate_kinematics(const double &res_dis_cut) {
  // common expression
  const double ME2 = 2 * effective_mass * neutrino.E();

  // determine max invariant mass (cannot be smaller than params::res_dis_cut)
  const double Wmax = min(res_dis_cut, sqrt(effective_mass2 + ME2) - lepton_mass);

  // choose random invariant mass (uniformly from [Wmin, Wmax])
  W = Wmin + (Wmax - Wmin) * frandom();
  W2 = W * W;

  // TODO: we integrate over z - what is its definition?
  const double z = frandom();

  // common expression
  const double W2_reduced = W2 - effective_mass2 - lepton_mass2;
  const double Mplus = effective_mass + 2 * neutrino.E();

  // aux variables
  const double A = (effective_mass + neutrino.E()) * W2_reduced + ME2 * neutrino.E();
  const double B = neutrino.E() * sqrt(pow2(W2_reduced - ME2) - 4 * lepton_mass2 * effective_mass * Mplus);
  const double C = 2 * effective_mass * Mplus;

  // energy transfer bounds
  const double q0_min = max((A - B) / C, lepton_mass);
  const double q0_max = min((A + B) / C, neutrino.E() - lepton_mass);

  // get random energy transfer
  q.t = q0_min + (q0_max - q0_min) * z * z;  // enhance low energy transfers are preferred

  // calculate jacobian
  jacobian = (q0_max - q0_min) * (Wmax - Wmin) * 2 * z;  // but compesated by this jakobian

  // temp kinematics variables
  const double q3 = sqrt(pow2(effective_mass + q.t) - W2);                                     // momentum transfer
  const double mom = sqrt(pow2(neutrino.E() - q.t) - lepton_mass2);                            // final lepton momentum
  const double cosine = (pow2(neutrino.E()) + pow2(mom) - pow2(q3)) / 2 / neutrino.E() / mom;  // scattering angle

  if (abs(cosine) > 1) return false;  // impossible kinematics
  if (pow2(neutrino.E() - q.t) - lepton_mass2 < 0) return false;

  vec mom_dir;                               // the unit vector in the direction of scattered lepton
  kinfinder(neutrino.p(), mom_dir, cosine);  // TODO: change this function!!!

  // final lepton kinematics done
  lepton = vect(neutrino.E() - q.t, mom_dir * mom);

  // update momentrum transfer
  q = neutrino.p() - lepton.p();

  // calculate the speed of hadronic system
  hadron_speed = q / sqrt(W2 + q3 * q3);

  return true;
}

bool res_kinematics::generate_kinematics(double _Q2, double _W)
{
  // fix invariant mass (for tests)
  W = _W;
  W2 = W * W;

  // fix Q2 (for tests)
  q.t = (W2 - effective_mass2 + _Q2)/2/effective_mass;

  // calculate jacobian
  jacobian = 1. / 2 / res_kinematics::avg_nucleon_mass;

  // temp kinematics variables
  const double q3 = sqrt(pow2(effective_mass + q.t) - W2);                                     // momentum transfer
  const double mom = sqrt(pow2(neutrino.E() - q.t) - lepton_mass2);                            // final lepton momentum
  const double cosine = (pow2(neutrino.E()) + pow2(mom) - pow2(q3)) / 2 / neutrino.E() / mom;  // scattering angle

  if (abs(cosine) > 1) return false;  // impossible kinematics
  if (pow2(neutrino.E() - q.t) - lepton_mass2 < 0) return false;

  vec mom_dir;                               // the unit vector in the direction of scattered lepton
  kinfinder(neutrino.p(), mom_dir, cosine);  // TODO: change this function!!!

  // final lepton kinematics done
  lepton = vect(neutrino.E() - q.t, mom_dir * mom);

  // update momentrum transfer
  q = neutrino.p() - lepton.p();

  // calculate the speed of hadronic system
  hadron_speed = q / sqrt(W2 + q3 * q3);

  return true;
}

void res_kinematics::set_kinematics(event &e) {
  // set kinematics necessary to calculate cross section
  W = e.W();
  W2 = W * W;
  q = e.res_q;
  neutrino = e.res_nu;
  jacobian = e.res_jacobian;
}

double get_binding_energy(const params &p, particle& target, nucleus &t) {
  switch (p.nucleus_target) {
    case 0:  // free nucleon
      return 0;
    case 1:  // (global) Fermi gas
      return p.nucleus_E_b;
    case 2:  // local Fermi gas
      return t.Ef(target) + p.kaskada_w;
    case 3:  // Bodek-Ritchie; temporary prescription taken from GFG
      return p.nucleus_E_b;
    case 4:  // effective spectral function
      return binen(target.p(), p.nucleus_p, p.nucleus_n);
    case 5:  // deuterium
      return deuter_binen(target.p());
    case 6:  // effective potential
        assert ( !"For a moment effective potential cannot be used for RES" );
      return p.nucleus_E_b;
    default:
      return 0;
  }
  return 0;
}
