#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>
#include <random>
#include <omp.h>
#include <nlohmann/json.hpp>
#include <icecream.hpp>
#include "Norm.hpp"
#include "PrivRepGame.hpp"

constexpr Reputation B = Reputation::B, G = Reputation::G;
constexpr Action C = Action::C, D = Action::D;

struct ParametersWithMu {
  size_t n_init;
  size_t n_steps;
  size_t N;
  double q;
  double mu_percept;
  double benefit;
  double beta;
  double mu;
  uint64_t seed;
  ParametersWithMu() : n_init(1e4), n_steps(1e4), N(30), q(0.9), mu_percept(0.05), benefit(5.0), beta(1.0), mu(0.01), seed(123456789) {};
  EvolPrivRepGame::SimulationParameters ToEvolParams() const {
    return EvolPrivRepGame::SimulationParameters(N, n_init, n_steps, q, mu_percept, seed);
  }

  NLOHMANN_DEFINE_TYPE_INTRUSIVE_WITH_DEFAULT(ParametersWithMu, n_init, n_steps, N, q, mu_percept, benefit, beta, mu, seed);
};

int main(int argc, char** argv) {
  using namespace nlohmann;

  std::vector<std::string> args;
  json j = json::object();
  bool swap_gb = false;
  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == "-j" && i + 1 < argc) {
      std::ifstream fin(argv[++i]);
      // check if file exists
      if (fin) {
        fin >> j;
        fin.close();
      } else {
        std::cerr << argv[i] << std::endl;
        std::istringstream iss(argv[i]);
        iss >> j;
      }
    }
    else if (std::string(argv[i]) == "-s" && i + 1 < argc) {
      swap_gb = std::string(argv[++i]) == "true";
    }
    else if (std::string(argv[i]) == "-h") {
      std::cerr << "Usage: " << argv[0] << " [options] norm" << std::endl;
      std::cerr << "Options:" << std::endl;
      std::cerr << "  -j <json>  : json string or file name" << std::endl;
      std::cerr << "  -s <bool>  : swap G and B" << std::endl;
      std::cerr << "  -h         : print this message" << std::endl;
      return 0;
    }
    else {
      args.emplace_back(argv[i]);
    }
  }

  ParametersWithMu params = j.get<ParametersWithMu>();
  if (args.size() != 1) {
    std::cerr << "[Error] invalid number of arguments" << std::endl;
    std::cerr << "Usage: " << argv[0] << " [options] norm" << std::endl;
    return 1;
  }
  Norm n = Norm::ParseNormString(args.at(0), swap_gb);

  std::cerr << "params: " << nlohmann::json(params).dump(2) << std::endl;
  std::cerr << "norm: " << n.Inspect() << std::endl;

  auto start = std::chrono::high_resolution_clock::now();

  auto evo_param = params.ToEvolParams();

  EvolPrivRepGameFiniteMutationRateAllCAllD evol(n, evo_param);
  auto result = evol.CalculateEquilibrium(params.benefit, params.beta, params.mu);
  IC(result.OverallCooperationLevel(), result.OverallAbundances() );

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cerr << "Elapsed time: " << elapsed.count() << " s\n";

  for (size_t nf = 0; nf <= result.N; nf++) {
    for (size_t nc = 0; nc <= result.N - nf; nc++) {
      std::cout << nf << ' ' << nc << ' ' << result.frequency[nf][nc] << ' ' << result.cooperation_level[nf][nc] << std::endl;
    }
  }

  return 0;
}