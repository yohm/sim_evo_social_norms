#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <chrono>
#include <random>
#include <nlohmann/json.hpp>
#include <icecream.hpp>
#include "Norm.hpp"


void PrintProgress(double progress) {
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
vd_t SolveByRungeKutta(std::function<void(const vd_t&,vd_t&)>& func, const vd_t& init, double dt, size_t n_iter) {
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
    // normalize ht
    double sum_inv = 1.0 / sum;
    for (int i = 0; i < N; i++) { ht[i] *= sum_inv; }
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

  // x_dot must have size N
  std::function<void(const vd_t&,vd_t&)> calc_x_dot = [&p_fix,&alpha,r_mut,N](const vd_t& x, vd_t& x_dot) {
    for (size_t i = 0; i < N; i++) {
      double dx = 0.0;
      for (size_t j = 0; j < N; j++) {
        if (i == j) continue;
        dx += (1.0 - r_mut) * x[i] * x[j] * alpha[j][i] - r_mut * x[i] * p_fix[i][j] / N + r_mut * x[j] * p_fix[i][j] / N;
      }
      x_dot[i] = dx;
    }
  };

  const size_t T_max = 500;
  vd_t x(N, 0.0);
  // measure elapsed time
  {
    std::ofstream fout("timeseries.dat");
    for (double& xi : x) { xi = 1.0 / N; }
    size_t dt = T_max / 500;
    for (size_t t = 0; t < T_max; t++) {
      x = SolveByRungeKutta(calc_x_dot, x, 0.1, 10);
      if (t % dt == 0) {
        fout << t << ' ';
        double pc = 0.0;
        for (size_t i = 0; i < N; i++) { pc += x[i] * self_coop_levels[i]; }
        fout << pc << ' ';
        for (double xi : x) { fout << xi << ' '; }
        fout << std::endl;
        PrintProgress(static_cast<double>(t) / T_max);
      }
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
