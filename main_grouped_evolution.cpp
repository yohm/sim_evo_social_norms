#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <thread>
#include <chrono>
#include <functional>
#include <nlohmann/json.hpp>
#include "Norm.hpp"
#include "PrivRepGame.hpp"


constexpr Reputation B = Reputation::B, G = Reputation::G;
constexpr Action C = Action::C, D = Action::D;

struct SimulationParams {
  size_t n_init;
  size_t n_steps;
  size_t N;
  double q;
  double mu_percept;
  double benefit;
  double beta;
  uint64_t seed;
  SimulationParams() : n_init(1e4), n_steps(1e4), N(30), q(0.9), mu_percept(0.05), benefit(5.0), beta(1.0), seed(123456789) {};

  NLOHMANN_DEFINE_TYPE_INTRUSIVE_WITH_DEFAULT(SimulationParams, n_init, n_steps, N, q, mu_percept, benefit, beta, seed);
};

double SelfCoopLevel(const Norm& norm, const SimulationParams& params) {
  PrivateRepGame prg({{norm, params.N}}, params.seed);
  prg.Update(params.n_init, params.q, params.mu_percept, false);
  prg.ResetCounts();
  prg.Update(params.n_steps, params.q, params.mu_percept, true);
  return prg.NormCooperationLevels()[0][0];
}

void CalculateFixationProbs(const SimulationParams& params, std::vector<std::vector<double>> & p_fix, std::vector<double> & self_coop_levels) {
  auto idx = [] (const Norm& n) { return n.IDwithoutR2(); };

  // loop over Norm
  constexpr size_t N_NORMS = 4096;
  assert(p_fix.size() == N_NORMS);
  assert(self_coop_levels.size() == N_NORMS);

  std::vector<Norm> unique_norms;
  std::vector<size_t> norm_index(N_NORMS, 0);
  for (int i = 0; i < N_NORMS; i++) {
    Norm norm = Norm::ConstructFromIDwithoutR2(i);

    if ( i < idx(norm.SwapGB()) ) {
      norm_index[i] = idx(norm.SwapGB());
    }
    else {
      norm_index[i] = i;
      unique_norms.push_back(norm);
    }
  }

  EvolPrivRepGame::SimulationParameters evoparams({params.n_init, params.n_steps, params.q, params.mu_percept, params.seed});

  #pragma omp parallel for schedule(static)
  for (size_t i = 0; i < unique_norms.size(); i++) {
    const Norm& n1 = unique_norms[i];
    double pc = SelfCoopLevel(n1, params);
    self_coop_levels[idx(n1)] = pc;
    p_fix[idx(n1)][idx(n1)] = 1.0 / params.N;
  }

  std::vector<std::array<size_t,2>> ij_pairs{};
  for (size_t i = 0; i < unique_norms.size(); i++) {
    for (size_t j = i+1; j < unique_norms.size(); j++) {
      ij_pairs.push_back({i, j});
    }
  }

  // loop over ij_pairs
  #pragma omp parallel for schedule(static)
  for (const auto& [i,j] : ij_pairs) {
    // std::cerr << "i, j: " << i << ", " << j << std::endl;
    const Norm& n1 = unique_norms[i];
    const Norm& n2 = unique_norms[j];
    EvolPrivRepGame evol(params.N, std::vector<Norm>({n1, n2}), evoparams);
    auto rhos = evol.FixationProbabilities(params.benefit, params.beta);
    p_fix[idx(n1)][idx(n2)] = rhos[0][1];
    p_fix[idx(n2)][idx(n1)] = rhos[1][0];
  }

  // calculate non-unique-norms
  for (size_t i = 0; i < N_NORMS; i++) {
    size_t ni = norm_index[i];
    self_coop_levels[i] = self_coop_levels[ni];
    for (size_t j = 0; j < N_NORMS; j++) {
      size_t nj = norm_index[j];
      p_fix[i][j] = p_fix[ni][nj];
    }
  }
}


int main(int argc, char* argv[]) {
  // run evolutionary simulation in group-structured population
  // strategy space: deterministic strategies without R2 (~ 2000 strategies)

  using namespace nlohmann;

  std::vector<std::string> args;
  json j = json::object();
  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == "-j" && i + 1 < argc) {
      std::ifstream fin(argv[++i]);
      // check if file exists
      if (fin) {
        fin >> j;
        fin.close();
      }
      else {
        std::istringstream iss(argv[i]);
        iss >> j;
      }
    }
    else {
      std::cerr << "unknown option: " << argv[i] << std::endl;
      return 1;
    }
  }
  SimulationParams params = j.get<SimulationParams>();

  // calculate intra-group fixation probabilities
  constexpr size_t N_NORMS = 4096;
  std::vector<std::vector<double>> p_fix(N_NORMS, std::vector<double>(N_NORMS, 0.0));
  std::vector<double> self_coop_levels(N_NORMS, 0.0);

  // measure elapsed time
  auto start = std::chrono::system_clock::now();

  CalculateFixationProbs(params, p_fix, self_coop_levels);
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cerr << "elapsed time: " << elapsed_seconds.count() << "s\n";

  // print fixation probabilities and cooperation levels
  std::ofstream fout("fixation_probs.dat");
  for (size_t i = 0; i < N_NORMS; i++) {
    fout << self_coop_levels[i] << " ";
    for (size_t j = 0; j < N_NORMS; j++) {
      fout << p_fix[i][j] << " ";
    }
    fout << std::endl;
  }

  // run evolutionary simulation in group-structured population
  // [TODO] implement me

  return 0;
}

