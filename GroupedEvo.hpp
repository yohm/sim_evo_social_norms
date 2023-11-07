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

  vd_t TimeEvolutionODE(double benefit, double sigma_out, double r_mut, size_t T_max, double dt, std::ostream& fout) {
    const size_t N = N_NORMS;

    std::vector<std::vector<double>> alpha(N, std::vector<double>(N, 0.0));
    // alpha[i][j] : flow from i to j
    //   = p_fix[i][j] * p_inter[i][j] - p_fix[j][i] * p_inter[j][i]
    for (size_t i = 0; i < N; i++) {
      for (size_t j = i+1; j < N; j++) {
        alpha[i][j] = CalcAlpha(i, j, benefit, sigma_out);
        alpha[j][i] = -alpha[i][j];
      }
    }

    std::vector<double> mu_out(N, 0.0);
    for (size_t i = 0; i < N; i++) {
      mu_out[i] = MutationOutFlow(i);
    }

    // x_dot must have size N
    std::function<void(const vd_t&,vd_t&)> calc_x_dot = [this,&alpha,&mu_out,r_mut,N](const vd_t& x, vd_t& x_dot) {
      for (size_t i = 0; i < N; i++) {
        double dx = 0.0;
        double mut_in = 0.0;
        for (size_t j = 0; j < N; j++) {
          if (i == j) continue;
          dx += (1.0 - r_mut) * x[i] * x[j] * alpha[j][i];  // - r_mut * x[i] * p_fix[i][j] / static_cast<double>(N-1) + r_mut * x[j] * p_fix[j][i] / static_cast<double>(N-1);
          mut_in += x[j] * p_fix[j][i];
        }
        x_dot[i] = dx + r_mut * ( mut_in - x[i] * mu_out[i] ) / static_cast<double>(N-1);
      }
    };

    vd_t x(N, 0.0);
    const size_t n_iter = 100;
    const auto N_max = static_cast<size_t>( std::round((double)T_max / (dt*n_iter)) );
    for (double& xi : x) { xi = 1.0 / static_cast<double>(N); }
    size_t n_interval = N_max / 500;
    if (n_interval == 0) { n_interval = 1; }
    for (size_t n = 0; n < N_max; n++) {
      x = SolveByRungeKutta(calc_x_dot, x, dt, n_iter, n_iter);
      if (n % n_interval == 0) {
        fout << n*dt*n_iter << ' ';
        double pc = 0.0;
        for (size_t i = 0; i < N; i++) { pc += x[i] * self_coop_levels[i]; }
        fout << pc << ' ';

        for (double xi : x) { fout << xi << ' '; }
        fout << std::endl;

        std::vector<double> x_dot(x.size(), 0.0);
        calc_x_dot(x, x_dot);
        double x_dot_max = 0.0;
        for (size_t i = 0; i < N; i++) {
          if (x_dot[i] > x_dot_max) { x_dot_max = x_dot[i]; }
        }
        fout << x_dot_max << ' ';

        PrintProgress(static_cast<double>(n) / T_max);
      }
    }
    return x;
  }

  static void PrintProgress(double progress) {
    static auto start = std::chrono::high_resolution_clock::now();
    int barWidth = 70;
    std::cerr << "\33[2K[";
    int pos = static_cast<int>(barWidth * progress);
    for (int i = 0; i < barWidth; ++i) {
      if (i < pos) std::cerr << "=";
      else if (i == pos) std::cerr << ">";
      else std::cerr << " ";
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cerr << "] " << int(progress * 100.0) << " %,  " << elapsed.count() << " sec\r";
    std::cerr.flush();
    if (progress >= 1.0) {
      std::cerr << std::endl;
      return;
    }
  }
};

#endif