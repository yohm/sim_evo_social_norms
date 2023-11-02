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
};

#endif