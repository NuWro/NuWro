#ifndef _hybrid_gateway_h_
#define _hybrid_gateway_h_
#include <array>
// Hybrid model input
int hybrid_ABCDE(double El_inc, double Q2, double W, double leptonmass, double nucleonmass, double *costheta_pi_in, int N, int *params, double (*strucfuncs)[5], double (*Inclusive)[5]);
std::array<double, 5> get_lepton_vec(double El_inc, double Q2, double W, double leptonmass, double nucleonmass, int *params);
#endif
