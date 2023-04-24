#include <iostream>
#include <fstream>
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

void CalculateFixationProbs(const SimulationParams& params, std::vector<std::vector<double>>& p_fix, std::vector<double>& self_coop_levels) {

  // loop over Norm
  auto for_each_unique_strategy = [] (std::function<void(const Norm&)> f) {
    for (int j = 0; j < 256; j++) {
      AssessmentRule R1 = AssessmentRule::MakeDeterministicRule(j);
      AssessmentRule R2 = AssessmentRule::KeepRecipient();

      for (int i = 0; i < 16; i++) {
        ActionRule P = ActionRule::MakeDeterministicRule(i);
        Norm norm(R1, R2, P);
        if ( norm.ID() < norm.SwapGB().ID() ) {
          continue;
        }
        if (P == ActionRule::ALLC() && norm != Norm::AllC()) {
          continue;  // it is equivalent to AllC
        }
        if (P == ActionRule::ALLD() && norm != Norm::AllD()) {
          continue;  // it is equivalent to AllD
        }
        f(norm);
      }
    }
  };

  EvolPrivRepGame::SimulationParameters evoparams({params.n_init, params.n_steps, params.q, params.mu_percept, params.seed});

  auto idx = [] (const Norm& n) { return (n.Rd.ID() << 4) + n.P.ID(); };

  for_each_unique_strategy([&] (const Norm& n1) {

    double pc = SelfCoopLevel(n1, params);
    self_coop_levels[idx(n1)] = pc;

    for_each_unique_strategy([&] (const Norm& n2) {
      if (n1.ID() >= n2.ID()) return;  // avoid double counting
      std::cerr << idx(n1) << " -> " << idx(n2) << std::endl;
      // calculate fixation probability
      EvolPrivRepGame evol(params.N, {n1, n2}, evoparams);
      auto rhos = evol.FixationProbabilities(params.benefit, params.beta);
      p_fix[idx(n1)][idx(n2)] = rhos[0][1];
      p_fix[idx(n2)][idx(n1)] = rhos[1][0];
    });
  });
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
  size_t num_norms = 4096;
  std::vector<std::vector<double>> p_fix(num_norms, std::vector<double>(params.N, 0.0));
  std::vector<double> self_coop_levels(num_norms, 0.0);
  CalculateFixationProbs(params, p_fix, self_coop_levels);

  // print fixation probabilities and cooperation levels
  std::ofstream fout("fixation_probs.dat");
  for (size_t i = 0; i < num_norms; i++) {
    fout << self_coop_levels[i] << " ";
    for (size_t j = 0; j < num_norms; j++) {
      fout << p_fix[i][j] << " ";
    }
    fout << std::endl;
  }

  // run evolutionary simulation in group-structured population
  // [TODO] implement me

  return 0;
}

