/*
This function calculates RES event
RES region is defined in the hadronic mass from 1080 to the parameter res_dis_cut introduced in params.txt

The following model of the #Sect is postulated

XSect (nu + N -> mu + N' + pi)/dW d omega = [XSect (nu + N -> mu + N' + pi)/dW d omega ]_Delta * alfadelta (W) +
                                            [XSect (nu + N -> mu + N' + pi)/dW d omega ]_DIS_SPP* SPP(W) * alfados(W) +
                                            [XSect (nu + N -> mu + N' + pi)/dW d omega ]_DIS_nonSPP

The main idea behind is that DIS contribution simulates non-resonant part

The functions alfa and beta are a priori arbitrary but they have to satisfy the conditions

  alfadis(res_dis_cut)=1, ....alfadelta(res_dis_cut=0,............alfadis (1080)=0

alfadis and alfadelta are defined in alfa.cc

In principle for each channel the functions should be different. They should be fitted to experimental data

Channels are labelled by 4 integers corresponding to: nu/anu   cc/nc   proton/neutron    pi+/pi0/pi-
*/

#include "resevent2.h"
#include <TMCParticle.h>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "LeptonMass.h"
#include "alfa.h"
#include "charge.h"
#include "delta.h"
#include "dis2res.h"
#include "dis_cr_sec.h"
#include "event1.h"
#include "fragmentation.h"
#include "grv94_bodek.h"
#include "jednostki.h"
#include "masses.h"
#include "parameters.h"
#include "params.h"
#include "pauli.h"
#include "pdg_name.h"
#include "res_xsec.h"
#include "singlepion.h"
#include "vect.h"
#include "nucleus.h"

extern "C" int pycomp_(const int *);

void resevent2(params &p, event &e, nucleus &t, bool cc) {
  e.weight = 0;              // if kinmetically forbidden
  e.flag.res_delta = false;  // pion from background by default
  e.res_angrew = 1.0;        // no xsec scaling by default

  res_kinematics kin(e, t);  // kinematics variables

  // check threshold for pion production (otherwise left e.weight = 0)
  if (not kin.is_above_threshold()) return;

  // generate random kinematics (return false in the case of impossible kinematics)
  if (not kin.generate_kinematics(p.res_dis_cut)) return;

  // save for reweighting
  e.res_q = kin.q;
  e.res_nu = kin.neutrino;
  e.res_jacobian = kin.jacobian;

  // save final lepton (kin.lepton is in target rest frame so boost it first)
  particle final_lepton = kin.lepton;
  final_lepton.boost(kin.target.v());
  final_lepton.pdg = kin.neutrino.pdg + cc * (1 - 2.0 * (kin.neutrino.pdg > 0));

  e.out.push_back(final_lepton);

  res_xsec xsec(kin, cc);  // cross section handler - initialize SPP indices and DIS constribution

  if (not kin.is_above_pythia_threshold() || xsec.is_no_dis()) {
    xsec.set_xsec_nopythia(kin, p);           // initialize xsec struct with SPP parameters
    e.weight = xsec.get_total(kin.jacobian);  // save total cross section

    // final state particles
    particle final_pion, final_nucleon;

    final_pion.set_pdg_and_mass(xsec.get_pion_pdg());
    final_nucleon.set_pdg_and_mass(nukleon_out_(kin.W, kin.neutrino.pdg, kin.target.pdg, final_pion.pdg, cc));

    // produces 4-momenta of final pair: nucleon + pion
    kin2part(kin.W, final_nucleon.pdg, final_pion.pdg, final_nucleon, final_pion);

    // TODO: apply rewangle corection

    // boost back to LAB frame
    final_nucleon.p4() = final_nucleon.boost(kin.hadron_speed);
    final_nucleon.p4() = final_nucleon.boost(kin.target.v());

    final_pion.p4() = final_pion.boost(kin.hadron_speed);
    final_pion.p4() = final_pion.boost(kin.target.v());

    // save final state hadrons
    e.out.push_back(final_pion);
    e.out.push_back(final_nucleon);
  } else  // W above pythia threshold and fromdis > 0
  {
    // the algorithm starts from the production of PYTHIA event
    TPythia6 *pythia71 = get_pythia();

    int nof_particles = 0;       // number of particles in the final state
    Pyjets_t *pythia_particles;  // pythia particles placeholder

    // force at least 5 particles in the final state
    // TODO: including initial particles?
    while (nof_particles < 5) {
      hadronization(kin.neutrino.E(), kin.W, kin.q.t, kin.lepton_mass, kin.neutrino.pdg, kin.target.pdg, cc);
      pythia_particles = pythia71->GetPyjets();
      nof_particles = pythia71->GetN();
    }

    /*
    There are three different possible outcomes:
      a) spp event nof_particles = 5
      b) more inelastic event nof_particles > 5, typically 7
      c) single kaon production; this causes technical complications because nof_particles = 5 also in this case
    */

    // single pion production -> 5 particles including a pion
    if (nof_particles == 5 and (PDG::pion(pythia_particles->K[1][3]) or PDG::pion(pythia_particles->K[1][4]))) {
      // K[1][3] or K[1][4] is a pion
      const int pion_pdg = PDG::pion(pythia_particles->K[1][3]) ? pythia_particles->K[1][3] : pythia_particles->K[1][4];
      // PDG to SPP code
      const int t = pdg2spp(pion_pdg);
      /*
        t:
          0 for pi+
          1 for pi0
          2 for pi-

        finalcharge:
          2 for proton + pi+
          1 for proton + pi0 or neutron + pi+
          0 for proton + pi- or neutron + pi0
         -1 for neutron + pi-

        finalcharge + t:
          1 -> neutron + pi- or neutron + pi0 or neutron + pi+
          2 -> proton + pi- or proton + pi+ or proton + pi0
      */
      const int nucleon_pdg = xsec.final_charge + t == 1 ? PDG::pdg_neutron : PDG::pdg_proton;

      // initialize xsec struct according to Pythia result
      xsec.set_xsec(kin, p, pion_pdg, nucleon_pdg, e.in[0].t);

      e.weight = xsec.get_total(kin.jacobian);  // save total cross sections

      // randomly decide if SPP comes from Delta or DIS
      if (xsec.get_dis_fraction() > frandom())  // SPP from DIS
        save_pythia_particles(e, pythia_particles, nof_particles, kin);
      else  // SPP from Delta
      {
        e.flag.res_delta = true;             // mark pion as comming from Delta
        particle final_nucleon, final_pion;  // final particles placeholders

        // boost leptons from nu-N CMS to the hadronic CMS
        kin.neutrino.boost(-kin.hadron_speed);
        kin.lepton.boost(-kin.hadron_speed);

        // produces 4-momenta of final pair: nucleon + pion with density matrix information
        // returns angrew (keep this value for reweighting)
        e.res_angrew = kin4part(kin.neutrino, kin.lepton, kin.W, nucleon_pdg, pion_pdg, final_nucleon, final_pion,
                                p.delta_angular);

        e.weight *= e.res_angrew;  // reweight according to angular correlation

        // boost back to nu-N CMS frame and then to LAB frame
        kin.neutrino.boost(kin.hadron_speed);
        kin.neutrino.boost(kin.target.v());

        kin.lepton.boost(kin.hadron_speed);
        kin.lepton.boost(kin.target.v());

        final_nucleon = final_nucleon.boost(kin.hadron_speed);
        final_nucleon = final_nucleon.boost(kin.target.v());

        final_pion = final_pion.boost(kin.hadron_speed);
        final_pion = final_pion.boost(kin.target.v());

        // set final hadrons PDG codes
        final_nucleon.pdg = nucleon_pdg;
        final_pion.pdg = pion_pdg;

        // save final hadrons in out vector
        e.out.push_back(final_pion);
        e.out.push_back(final_nucleon);
      }
    } else  // more inelastic final state or single kaon production
    {
      e.weight = xsec.get_total(kin.jacobian);
      save_pythia_particles(e, pythia_particles, nof_particles, kin);
    }

    delete pythia71;
  }

  // set all outgoing particles position to target nucleon position
  for (int j = 0; j < e.out.size(); j++) e.out[j].r = e.in[1].r;
}

TPythia6 *get_pythia() {
  TPythia6 *pythia71 = new TPythia6();

  // Setting Pythia parameters - done by Jaroslaw Nowak

  // stable pi0
  pythia71->SetMDCY(pycomp_(&pizero), 1, 0);

  pythia71->SetMSTU(20, 1);  // advirsory warning for unphysical flavour switch off
  pythia71->SetMSTU(23, 1);  // It sets counter of errors at 0
  pythia71->SetMSTU(26, 0);  // no warnings printed

  // PARJ(32)(D=1GeV)
  // is, with quark masses added, used to define the minimum allowable enrgy of a colour singlet parton system
  pythia71->SetPARJ(33, 0.1);

  // PARJ(33)-PARJ(34)(D=0.8GeV, 1.5GeV) are, with quark masses added, used to define the remaining energy below which
  // the fragmentation of a parton system is stopped and two final hadrons formed.
  pythia71->SetPARJ(34, 0.5);
  pythia71->SetPARJ(35, 1.0);

  // PARJ(36) (D=2.0GeV) represents the dependence of the mass of final quark pair for defining the stopping point of
  // the fragmentation. Strongly corlated with PARJ(33-35)
  pythia71->SetPARJ(37, 1.);

  // MSTJ(17) (D=2) number of attemps made to find two hadrons that have a combined mass below the cluster mass and thus
  // allow a cluster to decay rather than collaps
  pythia71->SetMSTJ(18, 3);

  return pythia71;
}

void save_pythia_particles(event &e, Pyjets_t *pythia_particles, const int nof_particles, const res_kinematics &kin) {
  // loop over Pythia particles
  for (int i = 0; i < nof_particles; i++) {
    // i-th Pythia particle converted to NuWro format
    particle p = get_pythia_particle(pythia_particles, i, kin);

    e.temp.push_back(p);  // all particles are stored in temp vector

    // only stable (ks == 1) particles are stored in out vector
    if (p.ks == 1) e.out.push_back(p);
  }
}

particle get_pythia_particle(Pyjets_t *pythia_particles, const int particle_id, res_kinematics kin) {
  particle p;

  // assign particle's four-momentum
  p.t = pythia_particles->P[3][particle_id] * GeV;
  p.x = pythia_particles->P[0][particle_id] * GeV;
  p.y = pythia_particles->P[1][particle_id] * GeV;
  p.z = pythia_particles->P[2][particle_id] * GeV;

  // rotate the particle produced by PYTHIA acoording to the direction of the momentum transfer
  // in PYTHIA it is assumed that this direction is the Z axis
  rotation(p, kin.q);

  p = p.boost(kin.hadron_speed);  // boost back to the nu-N CMS frame
  p = p.boost(kin.target.v());    // boost back to tha LAB frame

  p.ks = pythia_particles->K[0][particle_id];     // HEP particle status
  p.pdg = pythia_particles->K[1][particle_id];    // particle pdf
  p.orgin = pythia_particles->K[2][particle_id];  // HEP particle origin

  return p;
}
