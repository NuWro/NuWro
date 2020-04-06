/*
This function calculates RES events from the Hybrid model
*/

#include "resevent_hybrid.h"
#include "params.h"
#include "event1.h"
#include "dis/res_kinematics.h"
#include "particle.h"
#include "vect.h"
#include "dis/LeptonMass.h"
#include "hybrid_RES.h"
#include "hybrid_RES2.h"
#include "hybrid/hybrid_gateway.h"
#include "TMatrixD.h"
#include "TDecompLU.h"
#include "TVectorD.h"
#include "TH1D.h"


void resevent_hybrid(params &p, event &e, bool cc) {      // free nucleon only!
  e.weight = 0;           // if kinematically forbidden

  res_kinematics kin(e);  // kinematics variables

  // check threshold for pion production (otherwise left e.weight = 0)
  if (not kin.is_above_threshold()) return;

  // generate random kinematics (return false in the case of impossible kinematics)
  if (not kin.generate_kinematics(1500)) return;
  //if (not kin.generate_kinematics(1500, p.Q2, p.W)) return;

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
  //vec kierunek = -hybrid_dir_from_adler(0, 0, kin.neutrino, kin.lepton);

  // specify the params needed for ABCDE (note strange order!)
  int params[4];
  params[0] = 1;                                   // only CC for now
  params[2] = (1 - 2.0 * (kin.neutrino.pdg > 0));  // helicity
  // params[3] is the target nucleon, params[1] is the decay channel

  // cross section function
  double (*hybrid_xsec)(res_kinematics*, int*, vect) = hybrid_dsdQ2dW_tab;
  //double (*hybrid_xsec)(res_kinematics*, int*, vect) = hybrid_dsdQ2dWdcth_tab;
  //double (*hybrid_xsec)(res_kinematics*, int*, vect) = hybrid_dsdQ2dWdcth;
  //double (*hybrid_xsec)(res_kinematics*, int*, vect) = hybrid_dsdQ2dWdOm;

  switch (final_charge) { // calculate the cross section with only the CMS variables
    case  2:  // pi+ + proton (nu_11)
      {params[3] = 1; params[1] = 1;
       final_pion.set_pdg_and_mass( PDG::pdg_piP ); final_nucleon.set_pdg_and_mass( PDG::pdg_proton );
       kin1part(kin.W, final_nucleon.pdg, final_pion.pdg, final_nucleon, final_pion, kierunek);
       xsec_pip = hybrid_xsec(&kin, params, final_pion);}
      xsec_inclusive = xsec_pip;
      if ( not (xsec_inclusive > 0) ) return;
      break;
    case  1:  // pi+ + neutron (nu_22) or pi0 + proton (nu_21)
      {params[3] = 2; params[1] = 2;
       final_pion.set_pdg_and_mass( PDG::pdg_piP ); final_nucleon.set_pdg_and_mass( PDG::pdg_neutron );
       kin1part(kin.W, final_nucleon.pdg, final_pion.pdg, final_nucleon, final_pion, kierunek);
       xsec_pip = hybrid_xsec(&kin, params, final_pion);}
      {params[3] = 2; params[1] = 1;
       final_pion2.set_pdg_and_mass( PDG::pdg_pi ); final_nucleon2.set_pdg_and_mass( PDG::pdg_proton );
       kin1part(kin.W, final_nucleon2.pdg, final_pion2.pdg, final_nucleon2, final_pion2, kierunek);
       xsec_pi0 = hybrid_xsec(&kin, params, final_pion2);}
      xsec_inclusive = xsec_pip + xsec_pi0;
      if ( not (xsec_inclusive > 0) ) return;
      if( xsec_pip / xsec_inclusive < frandom() ) // random selection, switch to "2"
      {
        final_pion = final_pion2;
        final_nucleon = final_nucleon2;
        params[3] = 2; params[1] = 1;
      }
      else // make sure the params are okay for "1"
      {
        params[3] = 2; params[1] = 2;
      }
      break;
    case  0:  // pi0 + neutron (anu_11) or pi- + proton (anu_12)
      {params[3] = 1; params[1] = 1;
       final_pion.set_pdg_and_mass( PDG::pdg_pi ); final_nucleon.set_pdg_and_mass( PDG::pdg_neutron );
       kin1part(kin.W, final_nucleon.pdg, final_pion.pdg, final_nucleon, final_pion, kierunek);
       xsec_pi0 = hybrid_xsec(&kin, params, final_pion);}
      {params[3] = 1; params[1] = 2;
       final_pion.set_pdg_and_mass( -PDG::pdg_piP ); final_nucleon.set_pdg_and_mass( PDG::pdg_proton );
       kin1part(kin.W, final_nucleon2.pdg, final_pion2.pdg, final_nucleon2, final_pion2, kierunek);
       xsec_pim = hybrid_xsec(&kin, params, final_pion2);}
      xsec_inclusive = xsec_pi0 + xsec_pim;
      if ( not (xsec_inclusive > 0) ) return;
      if( xsec_pi0 / xsec_inclusive < frandom() ) // random selection, switch to "2"
      {
        final_pion = final_pion2;
        final_nucleon = final_nucleon2;
        params[3] = 1; params[1] = 2;
      }
      else // make sure the params are okay for "1"
      {
        params[3] = 1; params[1] = 1;
      }
      break;
    case -1:  // pi- + neutron (anu_22)
      {params[3] = 2; params[1] = 2;
       final_pion.set_pdg_and_mass( -PDG::pdg_piP ); final_nucleon.set_pdg_and_mass( PDG::pdg_neutron );
       kin1part(kin.W, final_nucleon.pdg, final_pion.pdg, final_nucleon, final_pion, kierunek);
       xsec_pim = hybrid_xsec(&kin, params, final_pion);}
      xsec_inclusive = xsec_pim;
      if ( not (xsec_inclusive > 0) ) return;
      break;
    default:
      cerr << "[WARNING]: Reaction charge out of range\n";
  };

  // Omega_pi^* was chosen in hadronic CMS

  // // Choose cos_theta^* for dsdQ2dW, in Adler frame
  // double costh_rnd = hybrid_sample_costh(kin.neutrino.E(), -kin.q*kin.q, kin.W, params);

  // // Choose phi^* for dsdQ2dWdcosth, in Adler frame
  // double phi_rnd = hybrid_sample_phi(kin.neutrino.E(), -kin.q*kin.q, kin.W, params, costh_rnd);

  // // Modify final hadron directions as specified in the Adler frame
  // vect nu  = kin.neutrino; nu.boost(-kin.hadron_speed); 
  // vect lep = kin.lepton;  lep.boost(-kin.hadron_speed);
  // vec dir_rnd = hybrid_dir_from_adler(costh_rnd, phi_rnd, nu, lep);

  // // Recalculate the hadronic kinematics, dir_rnd is the new direction of pion
  // double momentum = final_nucleon.momentum();
  // final_nucleon = vect(final_nucleon.E(), -momentum * dir_rnd.x, -momentum * dir_rnd.y, -momentum * dir_rnd.z);
  // final_pion = vect(final_pion.E(), momentum * dir_rnd.x, momentum * dir_rnd.y, momentum * dir_rnd.z);

  // set event weight
  e.weight = xsec_inclusive;

  // the cross section needs a factor (initial lepton momentum)^-2 in LAB
  // the cross section needs a jacobian: dw = dQ2/2M
  e.weight /= e.in[0].E() * e.in[0].E() / 2 / res_kinematics::avg_nucleon_mass / kin.jacobian;
  // coupling
  e.weight *= G*G*cos2thetac/2;
  // units
  e.weight /= cm2;

  // boost back to LAB frame
  final_nucleon.p4() = final_nucleon.boost(kin.hadron_speed);
  final_nucleon.p4() = final_nucleon.boost(kin.target.v());

  final_pion.p4() = final_pion.boost(kin.hadron_speed);
  final_pion.p4() = final_pion.boost(kin.target.v());

  // save final state hadrons
  // warning: the order is essential
  e.out.push_back(final_pion);
  e.out.push_back(final_nucleon);

  // set all outgoing particles position to target nucleon position
  for (int j = 0; j < e.out.size(); j++) e.out[j].r = e.in[1].r;
}

double hybrid_dsdQ2dW_tab(res_kinematics *kin, int params[4], vect final_pion)
{
  double result = 0.;
  double Q2 =-kin->q*kin->q;
  double W  = kin->W;

  double *hybrid_grid[] = {hybrid_grid_dQ2dW_11,   hybrid_grid_dQ2dW_22,   hybrid_grid_dQ2dW_21,
                           hybrid_grid_dQ2dW_2_11, hybrid_grid_dQ2dW_2_22, hybrid_grid_dQ2dW_2_21};
  int hg_idx = hybrid_grid_idx(kin->neutrino.pdg, params[3]*10 + params[1]);

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

  // we build leptonic tensor in CMS with q along z, see JES paper App. A
  vect kl_inc_lab = kin->neutrino.p4();   //
  vect kl_lab     = kin->lepton.p4();     //  They are not in LAB, but in target rest frame!
  vect q_lab      = kin->q;               //

  double v = q_lab.length() / (q_lab[0] + res_kinematics::avg_nucleon_mass);
  double g = 1 / sqrt(1 - v*v);
  double c = (kl_inc_lab.length() - kl_lab[3]) / q_lab.length();
  double s = sqrt(pow(kl_lab.length(),2) - pow(kl_lab[3],2)) / q_lab.length();

  vect kl_inc (g*(kl_inc_lab[0] - v*kl_inc_lab.length()*c), kl_inc_lab.length()*s,
               0, g*(-v*kl_inc_lab[0] + kl_inc_lab.length()*c));
  vect kl     (g*(kl_lab[0] - v*(kl_inc_lab.length()*c - q_lab.length())), kl_inc_lab.length()*s,
               0, g*(-v*kl_lab[0] + kl_inc_lab.length()*c - q_lab.length()));

  double kl_inc_dot_kl = kl_inc * kl;

  // cut for comparisons
  //if( Q2 > 1.91*GeV2 || W > 1400 ) return 0;
  //if( W > 1400 ) return 0;

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

    // 4 points surrounding the desired point
    int p00 = (Wf*Q2bin+Q2f)*5; // bottom left
    int p10 = p00+5;            // bottom right
    int p01 = p00+Q2bin*5;      // top left
    int p11 = p01+5;            // top right

    // interpolate
    for( int i = 0; i < 5; i++ )
      w[i] = bilinear_interp(hybrid_grid[hg_idx][p00+i], hybrid_grid[hg_idx][p10+i],
                             hybrid_grid[hg_idx][p01+i], hybrid_grid[hg_idx][p11+i], Q2d, Wd);
  }

  // calculate the leptonic tensor elements
  double l[5] = {0,0,0,0,0}; // 00, 03, 33, 11+22, 12
  l[0] = (2*kl_inc[0]*kl[0] - kl_inc_dot_kl);
  l[1] = ( -kl_inc[0]*kl[3] - kl[0]*kl_inc[3]);
  l[2] = (2*kl_inc[3]*kl[3] + kl_inc_dot_kl);
  l[3] = (2*kl_inc[1]*kl[1] + kl_inc_dot_kl);
  l[3]+= (2*kl_inc[2]*kl[2] + kl_inc_dot_kl);
  l[4] = (  kl_inc[0]*kl[3] - kl[0]*kl_inc[3]);

  // contract the tensors
  result = l[0]*w[0] + 2*l[1]*w[1] + l[2]*w[2] + 0.5*l[3]*w[3] - 2*l[4]*w[4];

  // correct the pion mass
  double W2 = W*W; double W4 = W2*W2;
  double nukmass2 = 881568.315821;
  double pionmass2 = 19054.761849;
  result /= sqrt(W4-2*W2*(pionmass2+nukmass2)+(nukmass2-pionmass2)*(nukmass2-pionmass2))/2.0/W;
  result *= final_pion.length();

  return result;
}

double hybrid_dsdQ2dWdcth_tab(res_kinematics *kin, int params[4], vect final_pion)
{
  double result = 0.;
  double Q2 =-kin->q*kin->q;
  double W  = kin->W;

  double *hybrid_grid[] = {hybrid_grid_dQ2dWdcth_11,   hybrid_grid_dQ2dWdcth_22,   hybrid_grid_dQ2dWdcth_21,
                           hybrid_grid_dQ2dWdcth_2_11, hybrid_grid_dQ2dWdcth_2_22, hybrid_grid_dQ2dWdcth_2_21};
  int hg_idx = hybrid_grid_idx(kin->neutrino.pdg, params[3]*10 + params[1]);

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
  double pion_momentum = final_pion.length();
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
  double pion_cos_theta = Zast * vec(final_pion) / pion_momentum;

  // we build leptonic tensor in CMS with q along z, see JES paper App. A
  vect kl_inc_lab = kin->neutrino.p4();   //
  vect kl_lab     = kin->lepton.p4();     //  They are not in LAB, but in target rest frame!
  vect q_lab      = kin->q;               //

  double v = q_lab.length() / (q_lab[0] + res_kinematics::avg_nucleon_mass);
  double g = 1 / sqrt(1 - v*v);
  double c = (kl_inc_lab.length() - kl_lab[3]) / q_lab.length();
  double s = sqrt(pow(kl_lab.length(),2) - pow(kl_lab[3],2)) / q_lab.length();

  vect kl_inc (g*(kl_inc_lab[0] - v*kl_inc_lab.length()*c), kl_inc_lab.length()*s,
               0, g*(-v*kl_inc_lab[0] + kl_inc_lab.length()*c));
  vect kl     (g*(kl_lab[0] - v*(kl_inc_lab.length()*c - q_lab.length())), kl_inc_lab.length()*s,
               0, g*(-v*kl_lab[0] + kl_inc_lab.length()*c - q_lab.length()));

  double kl_inc_dot_kl = kl_inc * kl;

  // cut for comparisons
  //if( Q2 > 1.91*GeV2 || W > 1400 ) return 0;

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

    // 8 points surrounding the desired point
    int p000 = (cthf*Wbin*Q2bin+Wf*Q2bin+Q2f)*5;       // bottom left close
    int p100 = p000+5;                                 // bottom right close
    int p010 = p000+Q2bin*5;                           // top left close
    int p110 = p010+5;                                 // top right close
    int p001 = p000+Wbin*Q2bin*5;                      // bottom left far
    int p101 = p001+5;                                 // bottom right far
    int p011 = p001+Q2bin*5;                           // top left far
    int p111 = p011+5;                                 // top right far

    // interpolate
    for( int i = 0; i < 5; i++ )
      w[i] = trilinear_interp(hybrid_grid[hg_idx][p000+i], hybrid_grid[hg_idx][p100+i],
                              hybrid_grid[hg_idx][p010+i], hybrid_grid[hg_idx][p110+i],
                              hybrid_grid[hg_idx][p001+i], hybrid_grid[hg_idx][p101+i],
                              hybrid_grid[hg_idx][p011+i], hybrid_grid[hg_idx][p111+i], Q2d, Wd, cthd);
  }

  // calculate the leptonic tensor elements
  double l[5] = {0,0,0,0,0}; // 00, 03, 33, 11+22, 12
  l[0] = (2*kl_inc[0]*kl[0] - kl_inc_dot_kl);
  l[1] = ( -kl_inc[0]*kl[3] - kl[0]*kl_inc[3]);
  l[2] = (2*kl_inc[3]*kl[3] + kl_inc_dot_kl);
  l[3] = (2*kl_inc[1]*kl[1] + kl_inc_dot_kl);
  l[3]+= (2*kl_inc[2]*kl[2] + kl_inc_dot_kl);
  l[4] = (  kl_inc[0]*kl[3] - kl[0]*kl_inc[3]);

  // contract the tensors
  result = l[0]*w[0] + 2*l[1]*w[1] + l[2]*w[2] + 0.5*l[3]*w[3] - 2*l[4]*w[4]; // *2Pi

  result *= pion_momentum / pow(2*Pi,3);                                      // /2Pi
  result *= 2; // Phase space!

  return result;
}

double hybrid_dsdQ2dWdcth(res_kinematics* kin, int params[4], vect final_pion)
{
  // placeholders
  double costh[1];
  double ABCDE[1][5] = {{0,0,0,0,0}};
  double result = 0.;

  // get Q^2, W
  double Q2 =-kin->q*kin->q;
  double W  = kin->W;

  // cut for comparisons
  //if( Q2 > 1.91*GeV2 || W > 1400 ) return 0;

  // find proper angles costh_pi^ast and phi_pi^ast
  double pion_momentum = final_pion.length();
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
  double pion_cos_theta = Zast * vec(final_pion) / pion_momentum;
  double pion_phi = atan2(Yast*vec(final_pion),Xast*vec(final_pion));

  // fill costh
  costh[0] = pion_cos_theta;

  // get ABCDE
  hybrid_ABCDE(kin->neutrino.E(), Q2, W, costh, 1, params, ABCDE);

  result  = ABCDE[0][0];                 // *2Pi
  result *= pion_momentum / pow(2*Pi,3); // /2Pi
  result *= 2; // Phase space!

  return result;
}

double hybrid_dsdQ2dWdOm(res_kinematics* kin, int params[4], vect final_pion)
{
  // placeholders
  double costh[1];
  double ABCDE[1][5] = {{0,0,0,0,0}};
  double result = 0.;

  // get Q^2, W
  double Q2 =-kin->q*kin->q;
  double W  = kin->W;

  // cut for comparisons
  //if( Q2 > 1.91*GeV2 || W > 1400 ) return 0;

  // find proper angles costh_pi^ast and phi_pi^ast
  double pion_momentum = final_pion.length();
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
  double pion_cos_theta = Zast * vec(final_pion) / pion_momentum;
  double pion_phi = atan2(Yast*vec(final_pion),Xast*vec(final_pion));

  // fill costh
  costh[0] = pion_cos_theta;

  // get ABCDE
  hybrid_ABCDE(kin->neutrino.E(), Q2, W, costh, 1, params, ABCDE);

  result = ABCDE[0][0] + ABCDE[0][1]*cos(pion_phi) + ABCDE[0][2]*cos(2*pion_phi)
                       + ABCDE[0][3]*sin(pion_phi) + ABCDE[0][4]*sin(2*pion_phi);
  result *= pion_momentum / pow(2*Pi,4);
  result *= 4 * Pi; // Phase space!

  return result;
}

double hybrid_sample_costh(double Enu, double Q2, double W, int params[4])
{
  // Choosing cos_th^* from dsdQ2dWdcosth
  double costh_rnd;

  // Specify points for the polynomial interpolation
  const int costh_pts = 3;             // 3, 5, 7, 9, ...
  double costh[costh_pts];             // Points of interpolation
  for( int i = 0; i < costh_pts; i++ ) // Fill costh with evenly spaced points
    costh[i] = -cos( Pi / (costh_pts-1) * i );

  // Get the A function, \propto dsdQ2dWdcosth
  double ABCDE[costh_pts][5]; double ds[costh_pts];
  hybrid_ABCDE(Enu, Q2, W, costh, costh_pts, params, ABCDE);
  for( int i = 0; i < costh_pts; i++ )
    ds[i] = ABCDE[i][0];

  // Fit a polynomial to given number of points in ds(costh)
  double poly_coeffs[costh_pts];
  hybrid_poly_fit(costh_pts, costh, ds, poly_coeffs);

  // Normalize the coefficients to obtain a probability density
  double norm = hybrid_poly_dist(costh_pts, poly_coeffs, costh[0], costh[costh_pts-1]);
  for( int i = 0; i < costh_pts; i++ )
    poly_coeffs[i] /= norm;

  // Choose a value for costh_rnd
  costh_rnd = hybrid_poly_rnd(costh_pts, poly_coeffs, costh[0], costh[costh_pts-1], 0.001);

  return costh_rnd;
}

double hybrid_sample_phi(double Enu, double Q2, double W, int params[4], double costh_rnd)
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
  hybrid_ABCDE(Enu, Q2, W, &costh_rnd, 1, params, ABCDE);
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

double hybrid_sample_phi_2(double Enu, double Q2, double W, int params[4], double costh_rnd)
{
  // Choosing phi^* from dsdQ2dWdOm
  double phi_rnd;

  // Get the ABCDE combination, \propto dsdQ2dWdOm
  double ABCDE[1][5];
  hybrid_ABCDE(Enu, Q2, W, &costh_rnd, 1, params, ABCDE);

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

double hybrid_sample_phi_3(double Enu, double Q2, double W, int params[4], double costh_rnd)
{
  // Choosing phi^* from dsdQ2dWdOm
  double phi_rnd;

  // Get the ABCDE combination, \propto dsdQ2dWdOm
  double ABCDE[1][5];
  hybrid_ABCDE(Enu, Q2, W, &costh_rnd, 1, params, ABCDE);

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

double hybrid_poly_rnd_2(const int N, double* coeffs, double x_min, double x_max, double epsilon)
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

  do
  {
    fx = ABCDE[0][0]*(x-x_min) + ABCDE[0][1]*(sin(x)-sin(x_min)) + ABCDE[0][2]*(sin(2*x)-sin(2*x_min))/2
                               - ABCDE[0][3]*(cos(x)-cos(x_min)) - ABCDE[0][4]*(cos(2*x)-cos(2*x_min))/2;
    fxp = ABCDE[0][0] + ABCDE[0][1]*cos(x) + ABCDE[0][2]*cos(2*x)
                      + ABCDE[0][3]*sin(x) + ABCDE[0][4]*sin(2*x);
    double new_x = x - (fx - y)/fxp;
    // if the function is very flat and the next step go out of bounds
    // use one bisective step to get out
    if( fabs(fxp) < epsilon || new_x < a || new_x > b )
    {
      if( fx > y )
        b = x;
      else
        a = x;
      x = (b-a)/2 + a;
    }
    else
    {
      x = new_x;
    }
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

void resevent_dir_hybrid(event& e)
{
  // get all 4-vectors from the event and boost them to the N-rest frame
  vect target   = e.in[1];  
  target.t     -= get_binding_energy (e.par, target);
  vect neutrino = e.in[0];  neutrino.boost(-target.v());
  vect lepton   = e.out[0]; lepton.boost  (-target.v());
  vect pion     = e.out[1]; pion.boost    (-target.v());
  vect nucleon  = e.out[2]; nucleon.boost (-target.v());

  // calculate hadron_speed and boost to CMS
  vect q = neutrino - lepton;
  double Mef = min(sqrt(target * target), res_kinematics::avg_nucleon_mass);
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
  double costh_rnd = hybrid_sample_costh(neutrino.t, -e.q2(), e.W(), params);

  // Choose phi^* for dsdQ2dWdcosth, in Adler frame. E_nu in N-rest!
  //double phi_rnd = hybrid_sample_phi(neutrino.t, -e.q2(), e.W(), params, costh_rnd);
  //double phi_rnd = hybrid_sample_phi_2(neutrino.t, -e.q2(), e.W(), params, costh_rnd);
  double phi_rnd = hybrid_sample_phi_3(neutrino.t, -e.q2(), e.W(), params, costh_rnd);

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

void resevent_phi_hybrid(event& e)
{
  // get all 4-vectors from the event and boost them to the N-rest frame
  vect target   = e.in[1];  
  target.t     -= get_binding_energy (e.par, target);
  vect neutrino = e.in[0];  neutrino.boost(-target.v());
  vect lepton   = e.out[0]; lepton.boost  (-target.v());
  vect pion     = e.out[1]; pion.boost    (-target.v());
  vect nucleon  = e.out[2]; nucleon.boost (-target.v());

  // calculate hadron_speed and boost to CMS
  vect q = neutrino - lepton;
  double Mef = min(sqrt(target * target), res_kinematics::avg_nucleon_mass);
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
  //double costh_rnd = hybrid_sample_costh(neutrino.t, -e.q2(), e.W(), params);
  vec Zast = neutrino - lepton;
  vec Yast = vecprod(neutrino,lepton);
  vec Xast = vecprod(Yast,Zast);
  Zast.normalize(); Yast.normalize(); Xast.normalize();
  double costh = Zast * vec(pion) / pion.length();

  // Choose phi^* for dsdQ2dWdcosth, in Adler frame. E_nu in N-rest!
  //double phi_rnd = hybrid_sample_phi(neutrino.t, -e.q2(), e.W(), params, costh_rnd);
  //double phi_rnd = hybrid_sample_phi_2(neutrino.t, -e.q2(), e.W(), params, costh_rnd);
  double phi_rnd = hybrid_sample_phi_3(neutrino.t, -e.q2(), e.W(), params, costh);

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

int hybrid_grid_idx(int neutrino_pdg, int channel)
{
  int hg_idx = -1;

  if(neutrino_pdg > 0)
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
      case 22: hg_idx = 0; break; // anu_22 has the tables of nu_11
      default:
        cerr << "[WARNING]: no hybrid inclusive grid!\n";
    }
  }

  return hg_idx;
}
