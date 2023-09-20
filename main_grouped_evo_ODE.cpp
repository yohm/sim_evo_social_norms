#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <chrono>
#include <random>
#include <nlohmann/json.hpp>
#include <icecream.hpp>
#include "Norm.hpp"
#include "EvolPrivRepGame.hpp"


nlohmann::json LoadMsgpackFile(const std::string& path) {
  std::ifstream fin(path, std::ios::binary);
  std::vector<char> bytes;
  // read the fin into the vector
  if (!fin) {
    std::cerr << "Error opening file " << path << std::endl;
    exit(1);
  }
  // get length of fin
  fin.seekg(0, fin.end);
  int length = fin.tellg();
  fin.seekg(0, fin.beg);

  // reserve space in vector for bytes
  bytes.resize(length);

  // read bytes into vector
  fin.read(bytes.data(), length);

    // unpack msgpack
  nlohmann::json j = nlohmann::json::from_msgpack(bytes);

  return j;
}

double InterGroupImitationProb(double pc_res, double pc_mut, double benefit, double sigma_out) {
  double pi_res = (benefit - 1.0) * pc_res;
  double pi_mut = (benefit - 1.0) * pc_mut;
  // f_{A\to B} = { 1 + \exp[ \sigma_out (s_A - s_B) ] }^{-1}
  return 1.0 / (1.0 + std::exp(sigma_out * (pi_res - pi_mut) ));
}

using vd_t = std::vector<double>;
vd_t SolveByRungeKutta(std::function<vd_t(vd_t)>& func, const vd_t& init, double dt, size_t n_iter) {
  const size_t N = init.size();
  vd_t ht = init;
  for (size_t t = 0; t < n_iter; t++) {
#ifdef DEBUG
    if (t % 10000 == 9999) {
      std::cerr << t << ' ' << ht[0] << ' ' << ht[1] << ' ' << ht[2] << std::endl;
    }
#endif
    vd_t k1 = func(ht);
    vd_t arg2(N, 0.0);
    for(int i = 0; i < k1.size(); i++) {
      k1[i] *= dt;
      arg2[i] = ht[i] + 0.5 * k1[i];
    }
    vd_t k2 = func(arg2);
    vd_t arg3(N, 0.0);
    for(int i = 0; i < k2.size(); i++) {
      k2[i] *= dt;
      arg3[i] = ht[i] + 0.5 * k2[i];
    }
    vd_t k3 = func(arg3);
    vd_t arg4(N, 0.0);
    for(int i = 0; i < k3.size(); i++) {
      k3[i] *= dt;
      arg4[i] = ht[i] + k3[i];
    }
    vd_t k4 = func(arg4);
    for(int i = 0; i < k4.size(); i++) {
      k4[i] *= dt;
    }
    vd_t delta(N, 0.0);
    double sum = 0.0;
    for (int i = 0; i < delta.size(); i++) {
      delta[i] = (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]) / 6.0;
      ht[i] += delta[i];
      sum += ht[i];
    }
    // normalize ht
    double sum_inv = 1.0 / sum;
    for (int i = 0; i < ht.size(); i++) { ht[i] *= sum_inv; }
  }
  return ht;
}

vd_t StationaryGroupedEvo(const std::vector<std::vector<double>>& p_fix, const std::vector<double>& self_coop_levels, double benefit, double sigma_out, double r_mut) {
  const size_t N = p_fix.size();

  std::vector<std::vector<double>> alpha(N, std::vector<double>(N, 0.0));
  // alpha[i][j] : flow from i to j
  //   = p_fix[i][j] * p_inter[i][j] - p_fix[j][i] * p_inter[j][i]
  for (size_t i = 0; i < N; i++) {
    for (size_t j = 0; j < N; j++) {
      double p_plus  = InterGroupImitationProb(self_coop_levels[i], self_coop_levels[j], benefit, sigma_out) * p_fix[i][j];
      double p_minus = InterGroupImitationProb(self_coop_levels[j], self_coop_levels[i], benefit, sigma_out) * p_fix[j][i];
      alpha[i][j] = p_plus - p_minus;
    }
  }

  std::function<vd_t(vd_t)> x_dot = [&p_fix,&alpha,r_mut,N](const vd_t& x) {
    vd_t ans(N, 0.0);
    for (size_t i = 0; i < N; i++) {
      double dx = 0.0;
      for (size_t j = 0; j < N; j++) {
        if (i == j) continue;
        dx += (1.0 - r_mut) * x[i] * x[j] * alpha[j][i] - r_mut * x[i] * p_fix[i][j] / N + r_mut * x[j] * p_fix[i][j] / N;
      }
      ans[i] = dx;
    }
    return ans;
  };

  const size_t T_max = 1;
  vd_t x(N, 0.0);
  // measure elapsed time
  {
    std::ofstream fout("timeseries.dat");
    for (double& xi : x) { xi = 1.0 / N; }
    size_t dt = T_max / 500;
    for (size_t t = 0; t < T_max; t++) {
      auto start = std::chrono::high_resolution_clock::now();
      x = SolveByRungeKutta(x_dot, x, 0.01, 100);
      if (t % dt == 0) {
        fout << t << ' ';
        double pc = 0.0;
        for (size_t i = 0; i < N; i++) { pc += x[i] * self_coop_levels[i]; }
        fout << pc << ' ';
        for (double xi : x) { fout << xi << ' '; }
        fout << std::endl;
      }
      auto end = std::chrono::high_resolution_clock::now();
      // get elapsed time in second
      std::chrono::duration<double> elapsed = end - start;

      std::cerr << "t: " << t << ", elapsed: " << elapsed.count() << std::endl;
    }
  }

  return x;
}

int main(int argc, char* argv[]) {
  // calculate stationary state by ODE for evolution in group-structured population

  using json = nlohmann::json;

  if (argc != 4) {
    std::cerr << "Usage: " << argv[0] << " <input_msgpack_file> <sigma_out> <relative mutation rate>" << std::endl;
    return 1;
  }

  // load fixation probabilities from the input file
  json j_in = LoadMsgpackFile(argv[1]);
  std::cerr << "json loaded" << std::endl;

  double sigma_out = std::stod(argv[2]);
  double relative_mutation_rate = std::stod(argv[3]);

  std::vector<Norm> norms;
  for (auto& norm_id_j : j_in["norm_ids"]) {
    int id = norm_id_j.get<int>();
    Norm n = Norm::ConstructFromID(id);
    norms.emplace_back(n);
  }
  std::cerr << "  # of norms: " << norms.size() << std::endl;

  const size_t N_NORMS = norms.size();
  std::vector<std::vector<double>> p_fix(N_NORMS, std::vector<double>(N_NORMS, 0.0));
  std::vector<double> self_coop_levels(N_NORMS, 0.0);
  for (size_t i = 0; i < N_NORMS; ++i) {
    for (size_t j = 0; j < N_NORMS; ++j) {
      p_fix[i][j] = j_in["p_fix"][i * N_NORMS + j].get<double>();
    }
  }
  for (size_t i = 0; i < self_coop_levels.size(); ++i) {
    self_coop_levels[i] = j_in["self_coop_levels"][i].get<double>();
  }
  double benefit = j_in["params"]["benefit"].get<double>();

  auto stationary = StationaryGroupedEvo(p_fix, self_coop_levels, benefit, sigma_out, relative_mutation_rate);

  // make tuple of norm_id & stationary & coop_level
  std::vector<std::tuple<int, double, double>> norms_to_measure;
  for (size_t i = 0; i < norms.size(); ++i) {
    norms_to_measure.emplace_back(norms[i].ID(), stationary[i], self_coop_levels[i]);
  }

  // sort by stationary
  std::sort(norms_to_measure.begin(), norms_to_measure.end(), [](auto& a, auto& b) {
    return std::get<1>(a) > std::get<1>(b);
  });

  // print out
  std::ofstream fout("stationary.dat");
  for (auto& n : norms_to_measure) {
    fout << std::get<0>(n) << ' ' << std::get<1>(n) << ' ' << std::get<2>(n) << std::endl;
  }

  return 0;
}
