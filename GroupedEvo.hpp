#ifndef GROUPED_EVO_HPP
#define GROUPED_EVO_HPP

#include <iostream>
#include <vector>

class GroupedEvo {
  public:
  GroupedEvo(const std::vector<size_t> &norm_ids, const std::vector<std::vector<double>> &p_fix,
             const std::vector<double> &self_coop_levels) :
    norm_ids(norm_ids), p_fix(p_fix), self_coop_levels(self_coop_levels) {
    N_NORMS = norm_ids.size();
    assert(p_fix.size() == N_NORMS);
    assert(self_coop_levels.size() == N_NORMS);
  }

  size_t N_NORMS;
  std::vector<size_t> norm_ids;
  std::vector<std::vector<double>> p_fix;
  std::vector<double> self_coop_levels;

  static double InterGroupImitationProb(double pc_res, double pc_mut, double benefit, double sigma_out) {
    double pi_res = (benefit - 1.0) * pc_res;
    double pi_mut = (benefit - 1.0) * pc_mut;
    // f_{A\to B} = { 1 + \exp[ \sigma_out (s_A - s_B) ] }^{-1}
    return 1.0 / (1.0 + std::exp(sigma_out * (pi_res - pi_mut)));
  }

  double CalcAlpha(size_t i, size_t j, double benefit, double sigma_out) const {
    // flow from i to j caused by inter-group imitation
    double p_plus = InterGroupImitationProb(self_coop_levels[i], self_coop_levels[j], benefit, sigma_out) * p_fix[i][j];
    double p_minus =
      InterGroupImitationProb(self_coop_levels[j], self_coop_levels[i], benefit, sigma_out) * p_fix[j][i];
    return p_plus - p_minus;
  }

  double MutationOutFlow(size_t i) {
    // sum_j p_fix[i][j]
    double mut_outflow = 0.0;
    for (size_t j = 0; j < N_NORMS; ++j) {
      mut_outflow += p_fix[i][j];
    }
    return mut_outflow;
  }

  using vd_t = std::vector<double>;
  static vd_t SolveByRungeKutta(std::function<void(const vd_t&,vd_t&)>& func, const vd_t& init, double dt, size_t n_iter, size_t normalize_interval = 0) {
    const size_t N = init.size();
    vd_t ht = init;
    vd_t k1(N, 0.0), arg2(N, 0.0), k2(N, 0.0), arg3(N, 0.0), k3(N, 0.0), arg4(N, 0.0), k4(N, 0.0);
    for (size_t t = 0; t < n_iter; t++) {
#ifdef DEBUG
      if (t % 10000 == 9999) {
      std::cerr << t << ' ' << ht[0] << ' ' << ht[1] << ' ' << ht[2] << std::endl;
    }
#endif
      func(ht, k1);
      for(int i = 0; i < N; i++) {
        k1[i] *= dt;
        arg2[i] = ht[i] + 0.5 * k1[i];
      }
      func(arg2, k2);
      for(int i = 0; i < N; i++) {
        k2[i] *= dt;
        arg3[i] = ht[i] + 0.5 * k2[i];
      }
      func(arg3, k3);
      for(int i = 0; i < N; i++) {
        k3[i] *= dt;
        arg4[i] = ht[i] + k3[i];
      }
      func(arg4, k4);
      for(int i = 0; i < N; i++) {
        k4[i] *= dt;
      }
      double sum = 0.0;
      for (int i = 0; i < N; i++) {
        ht[i] += (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]) / 6.0;
        sum += ht[i];
      }
      if (normalize_interval > 0 && t % normalize_interval == normalize_interval-1) {
        // normalize ht
        double sum_inv = 1.0 / sum;
        for (int i = 0; i < N; i++) { ht[i] *= sum_inv; }
      }
    }
    return ht;
  }
};

#endif