#include "hybrid/hybrid_gateway.h"
#include "mmapio.h"
#include "resevent_hybrid.h"
#include <array>
#include <cstring>

std::array<double, 5> get_nucleus_tensor(double W, double Q2, double costh,
                                         int *params) {
  double ABCDE[1][5], nuclear_tensor[1][5];
  hybrid_ABCDE(10 * GeV, Q2, W, 0.511 * MeV, res_kinematics::avg_nucleon_mass,
               &costh, 1, params, ABCDE, nuclear_tensor);
  return {nuclear_tensor[0][0] / 2, nuclear_tensor[0][1] / 2,
          nuclear_tensor[0][2] / 2, nuclear_tensor[0][3] / 2,
          nuclear_tensor[0][4] / 2};
};

int main(int argc, char **argv) {
  const size_t tab1_size = Q2bin_hybrid * Wbin_hybrid * cthbin_hybrid * 5;
  const size_t tab2_size = Q2bin_2_hybrid * Wbin_hybrid * cthbin_hybrid * 5;
  std::array<mmapio<double>, 3> tabulars{
      mmapio<double>{"hybrid_grid_dQ2dWdcth_11.dat", true, tab1_size},
      mmapio<double>{"hybrid_grid_dQ2dWdcth_22.dat", true, tab1_size},
      mmapio<double>{"hybrid_grid_dQ2dWdcth_21.dat", true, tab1_size},
  },
      tabular_dense{
          mmapio<double>{"hybrid_grid_dQ2dWdcth_2_11.dat", true, tab2_size},
          mmapio<double>{"hybrid_grid_dQ2dWdcth_2_22.dat", true, tab2_size},
          mmapio<double>{"hybrid_grid_dQ2dWdcth_2_21.dat", true, tab2_size}};
  // params [0]: channel (CC:1), [1]: decay channel, [2]: helicity (specify -1
  // here), [3]: target nucleon refer to hybrid_grid_idx
  // (src/resevent_hybrid.cc)
  int params_list[3][4] = {
      {1, 1, -1, 1} /*11*/, {1, 2, -1, 2} /*22*/, {1, 1, -1, 2} /*21*/};
  // enable parallelization if openmp founded,
// confirmed that parallelization is generating
// exact same result as serial version
#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for (size_t W_index = 0; W_index < Wbin_hybrid; W_index++) {
    double W = Wmin_hybrid + W_index * Wspc_hybrid;
    for (size_t channel_index = 0; channel_index < 3; channel_index++) {
      auto params = params_list[channel_index];
      for (size_t cth_index = 0; cth_index < cthbin_hybrid; cth_index++) {
        double cth = cthmin_hybrid + cth_index * cthspc_hybrid;
        for (size_t Q2_index = 0; Q2_index < Q2bin_hybrid; Q2_index++) {
          double Q = Q2min_hybrid + Q2_index * Q2spc_hybrid;
          auto nucleus_tensor = get_nucleus_tensor(W, Q, cth, params);
          // __builtin_trap();
          size_t index =
              index_calculator(cth_index, W_index, Q2_index, Q2bin_hybrid) * 5;
          memcpy(&tabulars[channel_index][index], nucleus_tensor.data(),
                 5 * sizeof(double));
        }
        for (size_t Q2_index = 0; Q2_index < Q2bin_2_hybrid; Q2_index++) {
          double Q = Q2min_2_hybrid + Q2_index * Q2spc_2_hybrid;
          auto nucleus_tensor = get_nucleus_tensor(W, Q, cth, params);
          size_t index =
              index_calculator(cth_index, W_index, Q2_index, Q2bin_2_hybrid) *
              5;
          memcpy(&tabular_dense[channel_index][index], nucleus_tensor.data(),
                 5 * sizeof(double));
        }
      }
    }
  }
  return 0;
}