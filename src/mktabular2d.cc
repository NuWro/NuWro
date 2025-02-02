#include "hybrid/hybrid_gateway.h"
#include "mmapio.h"
#include "resevent_hybrid.h"
#include <array>
#include <cstring>
#include <functional>

class nuclear_tensor_t : public std::array<double, 5> {
public:
  template <typename... T>
  nuclear_tensor_t(T &&...t) : std::array<double, 5>(std::forward<T>(t)...) {}
  nuclear_tensor_t() = default;
  nuclear_tensor_t(const nuclear_tensor_t &other) = default;
  nuclear_tensor_t(nuclear_tensor_t &&other) = default;
  nuclear_tensor_t operator/(double x) const {
    nuclear_tensor_t ret;
    for (size_t i = 0; i < 5; i++) {
      ret[i] = (*this)[i] / x;
    }
    return ret;
  }
  nuclear_tensor_t operator+(const nuclear_tensor_t &x) const {
    nuclear_tensor_t ret;
    for (size_t i = 0; i < 5; i++) {
      ret[i] = (*this)[i] + x[i];
    }
    return ret;
  }
  nuclear_tensor_t operator-(const nuclear_tensor_t &x) const {
    nuclear_tensor_t ret;
    for (size_t i = 0; i < 5; i++) {
      ret[i] = (*this)[i] - x[i];
    }
    return ret;
  }
  nuclear_tensor_t &operator=(const nuclear_tensor_t &other) = default;
  nuclear_tensor_t &operator=(nuclear_tensor_t &&other) = default;
  nuclear_tensor_t &operator+=(const nuclear_tensor_t &x) {
    for (size_t i = 0; i < 5; i++) {
      (*this)[i] += x[i];
    }
    return *this;
  }
  nuclear_tensor_t &operator-=(const nuclear_tensor_t &x) {
    for (size_t i = 0; i < 5; i++) {
      (*this)[i] -= x[i];
    }
    return *this;
  }
  nuclear_tensor_t &operator*=(double x) {
    for (size_t i = 0; i < 5; i++) {
      (*this)[i] *= x;
    }
    return *this;
  }
};

nuclear_tensor_t operator*(double a, const nuclear_tensor_t &b) {
  nuclear_tensor_t ret;
  for (size_t i = 0; i < 5; i++) {
    ret[i] = a * b[i];
  }
  return ret;
}

nuclear_tensor_t operator*(const nuclear_tensor_t &b, double a) {
  nuclear_tensor_t ret;
  for (size_t i = 0; i < 5; i++) {
    ret[i] = a * b[i];
  }
  return ret;
}

nuclear_tensor_t get_nucleus_tensor_3d(double W, double Q2, double costh,
                                       int *params) {
  double ABCDE[1][5], nuclear_tensor[1][5];
  hybrid_ABCDE(10 * GeV, Q2, W, 0.511 * MeV, res_kinematics::avg_nucleon_mass,
               &costh, 1, params, ABCDE, nuclear_tensor);
  return std::array<double, 5>{
      nuclear_tensor[0][0] / 2, nuclear_tensor[0][1] / 2,
      nuclear_tensor[0][2] / 2, nuclear_tensor[0][3] / 2,
      nuclear_tensor[0][4] / 2};
};

template <typename T>
T boole_integrate(std::function<T(double)> &func, double a, double b) {
  double h = (b - a) / 4;
  return (7 * func(a) + 32 * func(a + h) + 12 * func(a + 2 * h) +
          32 * func(a + 3 * h) + 7 * func(a + 4 * h)) *
         2 * h / 45;
}

template <typename T>
T simpson_integrate(std::function<T(double)> &func, double a, double b) {
  double h = (b - a) / 2;
  return h / 3 * (func(a) + 4 * func(a + h) + func(b));
}

template <typename T>
T do_integrate(std::function<T(double)> func, double a, double b, int divide) {
  double h = (b - a) / divide;
  T ret{};
  for (int i = 0; i < divide; i++) {
    ret += boole_integrate(func, a + i * h, a + (i + 1) * h);
  }
  return ret;
}

nuclear_tensor_t get_nucleus_tensor(double W, double Q2, int *params) {
  return do_integrate<nuclear_tensor_t>(
      [=](double costh) -> nuclear_tensor_t {
        return get_nucleus_tensor_3d(W, Q2, costh, params);
      },
      -1, 1, 128);
}

double get_pion_momentum(double hama) {
  double nukmass, pionmass;
  double W2 = hama * hama;
  double W4 = W2 * W2;

  // put average masses and return the momentum in that case
  double nukmass2 = 881568.315821;
  double pionmass2 = 19054.761849;
  return sqrt(W4 - 2 * W2 * (pionmass2 + nukmass2) +
              (nukmass2 - pionmass2) * (nukmass2 - pionmass2)) /
         2.0 / hama;
}

int main(int argc, char **argv) {
  const size_t tab1_size = Q2bin_hybrid * Wbin_hybrid * 5;
  const size_t tab2_size = Q2bin_2_hybrid * Wbin_hybrid * 5;
  std::array<mmapio<double>, 3> tabulars{
      mmapio<double>{"hybrid_grid_dQ2dW_11.dat", true, tab1_size},
      mmapio<double>{"hybrid_grid_dQ2dW_22.dat", true, tab1_size},
      mmapio<double>{"hybrid_grid_dQ2dW_21.dat", true, tab1_size},
  },
      tabular_dense{
          mmapio<double>{"hybrid_grid_dQ2dW_2_11.dat", true, tab2_size},
          mmapio<double>{"hybrid_grid_dQ2dW_2_22.dat", true, tab2_size},
          mmapio<double>{"hybrid_grid_dQ2dW_2_21.dat", true, tab2_size}};
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
      double pion_momentum = get_pion_momentum(W);
      for (size_t Q2_index = 0; Q2_index < Q2bin_hybrid; Q2_index++) {
        double Q = Q2min_hybrid + Q2_index * Q2spc_hybrid;
        auto nucleus_tensor = get_nucleus_tensor(W, Q, params);
        nucleus_tensor *= pion_momentum / pow(2 * Pi, 3);
        size_t index = index_calculator_2d(W_index, Q2_index, Q2bin_hybrid) * 5;
        memcpy(&tabulars[channel_index][index], nucleus_tensor.data(),
               5 * sizeof(double));
        // for (size_t g{}; g < 5; g++) {
        //   if (isnan(tabulars[channel_index][index + g]) ||
        //       isinf(tabulars[channel_index][index + g])) {
        //     std::cerr << "W: " << W << " Q2: " << Q
        //               << " channel: " << channel_index << " g: " << g
        //               << " value: " << tabulars[channel_index][index + g]
        //               << std::endl;
        //   }
        // }
      }
      for (size_t Q2_index = 0; Q2_index < Q2bin_2_hybrid; Q2_index++) {
        double Q = Q2min_2_hybrid + Q2_index * Q2spc_2_hybrid;
        auto nucleus_tensor = get_nucleus_tensor(W, Q, params);
        nucleus_tensor *= pion_momentum / pow(2 * Pi, 3);
        size_t index = index_calculator_2d(W_index, Q2_index, Q2bin_2_hybrid) * 5;
        memcpy(&tabular_dense[channel_index][index], nucleus_tensor.data(),
               5 * sizeof(double));
        // for (size_t g{}; g < 5; g++) {
        //   if (isnan(tabular_dense[channel_index][index + g]) ||
        //       isinf(tabular_dense[channel_index][index + g])) {
        //     std::cerr << "W: " << W << " Q2: " << Q
        //               << " channel: " << channel_index << " g: " << g
        //               << " value: " << tabular_dense[channel_index][index + g]
        //               << std::endl;
        //   }
        // }
      }
    }
  }
  return 0;
}