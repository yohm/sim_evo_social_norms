#include <iostream>
#include <fstream>
#include <cassert>
#include <chrono>
#include <regex>
#include <icecream.hpp>
#include <nlohmann/json.hpp>
#include "PrivRepGame.hpp"


template <typename T>
bool IsAllClose(T a, T b, double epsilon = 0.02) {
  for (size_t i = 0; i < a.size(); i++) {
    if (std::abs(a[i] - b[i]) > epsilon) {
      return false;
    }
  }
  return true;
}

bool IsClose(double a, double b, double epsilon = 0.02) {
  return std::abs(a - b) < epsilon;
}

void test_SelfCooperationLevel(const Norm& norm, double expected_c_level, double expected_good_rep) {
  auto start = std::chrono::high_resolution_clock::now();

  PrivateRepGame priv_game( {{norm, 50}}, 123456789ull);
  priv_game.Update(1e4, 0.9, 0.05, false);
  priv_game.ResetCounts();
  priv_game.Update(1e4, 0.9, 0.05, true);
  IC( priv_game.NormCooperationLevels(), priv_game.NormAverageReputation() );
  assert( IsClose(priv_game.SystemWideCooperationLevel(), expected_c_level, 0.02) );
  assert( IsAllClose(priv_game.NormCooperationLevels()[0], {expected_c_level}, 0.02) );
  assert( IsAllClose(priv_game.NormAverageReputation()[0], {expected_good_rep}, 0.02) );

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cerr << "Elapsed time: " << elapsed.count() << " s\n";
}

void test_RandomNorm() {
  // random norm
  test_SelfCooperationLevel(Norm::Random(), 0.5, 0.5);
  std::cerr << "test_RandomNorm passed" << std::endl;
}


void test_LeadingEight() {

  // measure time
  test_SelfCooperationLevel(Norm::L1(), 0.90, 0.90);
  std::cerr << "test L1 passed" << std::endl;

  test_SelfCooperationLevel(Norm::L2(), 0.66, 0.65);
  std::cerr << "test L2 passed" << std::endl;

  test_SelfCooperationLevel(Norm::L3(), 0.90, 0.90);
  std::cerr << "test L3 passed" << std::endl;

  test_SelfCooperationLevel(Norm::L4(), 0.90, 0.90);
  std::cerr << "test L4 passed" << std::endl;

  test_SelfCooperationLevel(Norm::L5(), 0.70, 0.70);
  std::cerr << "test L5 passed" << std::endl;

  test_SelfCooperationLevel(Norm::L6(), 0.50, 0.50);
  std::cerr << "test L6 passed" << std::endl;

  test_SelfCooperationLevel(Norm::L7(), 0.88, 0.88);
  std::cerr << "test L7 passed" << std::endl;

  test_SelfCooperationLevel(Norm::L8(), 0.0, 0.0);
  std::cerr << "test L8 passed" << std::endl;
}

void test_SelectionMutationEquilibrium() {
  auto start = std::chrono::high_resolution_clock::now();

  Norm norm = Norm::L1();
  EvolPrivRepGame::SimulationParameters params;
  params.n_init = 1e5;
  params.n_steps = 1e5;

  EvolPrivRepGame evol(50, {norm, Norm::AllC(), Norm::AllD()}, params);
  auto rhos = evol.FixationProbabilities(5.0, 1.0);
  IC(rhos);
  assert( IsClose(rhos[0][1], 0.097, 0.02) );
  assert( IsClose(rhos[0][2], 0.000, 0.02) );
  assert( IsClose(rhos[1][0], 0.012, 0.02) );
  assert( IsClose(rhos[2][0], 0.043, 0.02) );
  auto eq = evol.EquilibriumPopulationLowMut(rhos);
  IC(eq);
  assert(IsAllClose(eq, {0.30, 0.04, 0.66}, 0.02) );

  std::cerr << "test_SelectionMutationEquilibrium passed" << std::endl;

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cerr << "Elapsed time: " << elapsed.count() << " s\n";
}

void test_SelectionMutationEquilibrium2() {
  auto start = std::chrono::high_resolution_clock::now();

  EvolPrivRepGame::SimulationParameters params;
  params.n_init = 1e5;
  params.n_steps = 1e5;

  EvolPrivRepGameAllCAllD evol(50, params, 5.0, 1.0);

  auto selfc_rho_eq = evol.EquilibriumCoopLevelAllCAllD(Norm::L1());
  double self_cooperation_level = std::get<0>(selfc_rho_eq);
  auto rhos = std::get<1>(selfc_rho_eq);
  auto eq = std::get<2>(selfc_rho_eq);

  IC(self_cooperation_level);
  assert( IsClose(self_cooperation_level, 0.90, 0.02) );
  IC(rhos);
  assert( IsClose(rhos[0][1], 0.097, 0.02) );
  assert( IsClose(rhos[0][2], 0.000, 0.02) );
  assert( IsClose(rhos[1][0], 0.012, 0.02) );
  assert( IsClose(rhos[2][0], 0.043, 0.02) );
  IC(eq);
  assert(IsAllClose(eq, {0.30, 0.04, 0.66}, 0.02) );

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cerr << "Elapsed time: " << elapsed.count() << " s\n";
}

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


void CompareWithLocalMutants(const Norm& norm, const SimulationParams& params);

void PrintSelectionMutationEquilibriumAllCAllD(const Norm& norm, const SimulationParams& params, bool check_local_mutants = false) {
  auto start = std::chrono::high_resolution_clock::now();

  PrivateRepGame prg({{norm, params.N}}, params.seed);
  prg.Update(params.n_init, params.q, params.mu_percept, false);
  prg.ResetCounts();
  prg.Update(params.n_steps, params.q, params.mu_percept, true);
  IC( prg.NormAverageReputation(), prg.NormCooperationLevels());

  EvolPrivRepGame::SimulationParameters evo_params(params.n_init, params.n_steps, params.q, params.mu_percept, params.seed);
  EvolPrivRepGameAllCAllD evol(params.N, evo_params, params.benefit, params.beta);

  auto selfc_rho_eq = evol.EquilibriumCoopLevelAllCAllD(norm);
  double self_cooperation_level = std::get<0>(selfc_rho_eq);
  auto rhos = std::get<1>(selfc_rho_eq);
  auto eq = std::get<2>(selfc_rho_eq);

  IC(self_cooperation_level, rhos, eq);

  if (check_local_mutants) {
    CompareWithLocalMutants(norm, params);
  }

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cerr << "Elapsed time: " << elapsed.count() << " s\n";
}

void PrintCompetition(const Norm& n1, const Norm& n2, const SimulationParams& params) {
  auto start = std::chrono::high_resolution_clock::now();

  EvolPrivRepGame::SimulationParameters evo_params(params.n_init, params.n_steps, params.q, params.mu_percept, params.seed);

  EvolPrivRepGame evol(params.N, {n1, n2}, evo_params);
  auto fixs = evol.FixationProbabilities(params.benefit, params.beta);
  auto eq = evol.EquilibriumPopulationLowMut(fixs);
  IC( fixs, eq );

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cerr << "Elapsed time: " << elapsed.count() << " s\n";
}

void CompareWithLocalMutants(const Norm& norm, const SimulationParams& params) {
  EvolPrivRepGame::SimulationParameters evo_params(params.n_init, params.n_steps, params.q, params.mu_percept, params.seed);

  for (int i = 0; i < 20; i++) {
    auto serialized = norm.Serialize();
    const double delta = 1.0;
    if (serialized[i] <= 0.5) {
      serialized[i] += delta;
      if (serialized[i] > 1.0) {
        serialized[i] = 1.0;
      }
    }
    else {
      serialized[i] -= delta;
      if (serialized[i] < 0.0) {
        serialized[i] = 0.0;
      }
    }
    Norm mutant = Norm::FromSerialized(serialized);
    std::cout << mutant.Inspect();
    EvolPrivRepGame evol(params.N, {norm, mutant}, evo_params);
    auto rhos = evol.FixationProbabilities(params.benefit, params.beta);
    auto eq = evol.EquilibriumPopulationLowMut(rhos);
    IC(eq);
  }
}

Norm ParseNorm(const std::string& str) {
  std::regex re_d(R"(\d+)"); // regex for digits
  std::regex re_x(R"(^0x[0-9a-fA-F]+$)");  // regex for digits in hexadecimal
  // regular expression for 20 floating point numbers separated by space
  std::regex re_a(R"(^(0|1)(\.\d+)?( (0|1)(\.\d+)?){19}$)");
  if (std::regex_match(str, re_d)) {
    int id = std::stoi(str);
    return Norm::ConstructFromID(id);
  }
  else if (std::regex_match(str, re_x)) {
    int id = std::stoi(str, nullptr, 16);
    return Norm::ConstructFromID(id);
  }
  else if (std::regex_match(str, re_a)) {
    std::istringstream iss(str);
    std::array<double,20> serialized = {};
    for (int i = 0; i < 20; ++i) {
      iss >> serialized[i];
    }
    return Norm::FromSerialized(serialized);
  }
  else {
    return Norm::ConstructFromName(str);
  }
}

int main(int argc, char *argv[]) {

  std::vector<std::string> args;
  bool check_local_mutants = false;
  nlohmann::json j = nlohmann::json::object();
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
    else if (std::string(argv[i]) == "-l") {
      check_local_mutants = true;
    }
    else {
      args.push_back(argv[i]);
    }
  }

  if (args.size() == 0) {
    test_RandomNorm();
    test_LeadingEight();
    test_SelectionMutationEquilibrium();
    test_SelectionMutationEquilibrium2();
  }
  else if (args.size() == 1) {
    Norm n = ParseNorm(args.at(0));
    std::cout << n.Inspect();
    SimulationParams params = j.get<SimulationParams>();
    std::cout << nlohmann::json(params).dump(2) << std::endl;
    PrintSelectionMutationEquilibriumAllCAllD(n, params, check_local_mutants);
  }
  else if (args.size() == 2) {  // if two arguments are given, direct competition between two norms are shown
    Norm n1 = ParseNorm(args.at(0));
    Norm n2 = ParseNorm(args.at(1));
    SimulationParams params = j.get<SimulationParams>();
    std::cout << nlohmann::json(params).dump(2) << std::endl;
    PrintCompetition(n1, n2, params);
  }

  return 0;
}