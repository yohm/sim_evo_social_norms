#include <iostream>
#include <fstream>
#include <cassert>
#include <chrono>
#include <regex>
#include <icecream.hpp>
#include <nlohmann/json.hpp>
#include "PrivRepGame.hpp"

#include <gtest/gtest.h>

void test_SelfCooperationLevel(const Norm& norm, double expected_c_level, double expected_good_rep) {
  PrivateRepGame priv_game( {{norm, 50}}, 123456789ull);
  priv_game.Update(1e4, 0.9, 0.05, false);
  priv_game.ResetCounts();
  priv_game.Update(1e4, 0.9, 0.05, true);
  // IC( priv_game.NormCooperationLevels(), priv_game.NormAverageReputation() );
  EXPECT_NEAR( priv_game.SystemWideCooperationLevel(), expected_c_level, 0.02);
  EXPECT_NEAR( priv_game.NormCooperationLevels()[0][0], expected_c_level, 0.02);
  EXPECT_NEAR( priv_game.NormAverageReputation()[0][0], expected_good_rep, 0.02);
}

TEST(SelfCooperationLevel, RandomNorm) {
  test_SelfCooperationLevel(Norm::Random(), 0.5, 0.5);
}

TEST(SelfCooperationLevel, LeadingEight) {
  test_SelfCooperationLevel(Norm::L1(), 0.90, 0.90);
  test_SelfCooperationLevel(Norm::L2(), 0.66, 0.65);
  test_SelfCooperationLevel(Norm::L3(), 0.90, 0.90);
  test_SelfCooperationLevel(Norm::L4(), 0.90, 0.90);
  test_SelfCooperationLevel(Norm::L5(), 0.70, 0.70);
  test_SelfCooperationLevel(Norm::L6(), 0.50, 0.50);
  test_SelfCooperationLevel(Norm::L7(), 0.88, 0.88);
  test_SelfCooperationLevel(Norm::L8(), 0.0, 0.0);
}

TEST(SelectionMutationEquilibrium, L1_AllC_AllD) {
  Norm norm = Norm::L1();
  EvolPrivRepGame::SimulationParameters params;
  params.n_init = 1e5;
  params.n_steps = 1e5;

  EvolPrivRepGame evol(50, {norm, Norm::AllC(), Norm::AllD()}, params);
  auto rhos = evol.FixationProbabilities(5.0, 1.0);
  // rhos: [
  // [0, 0.0976933, 4.47544e-29],
  // [0.0114276, 0, 0.668377],
  // [0.0444453, 2.16425e-24, 0]
  // ]

  EXPECT_NEAR(rhos[0][1], 0.097, 0.02);
  EXPECT_NEAR(rhos[0][2], 0.000, 0.02);
  EXPECT_NEAR(rhos[1][0], 0.012, 0.02);
  EXPECT_NEAR(rhos[2][0], 0.043, 0.02);
  EXPECT_NEAR(rhos[1][2], 0.668, 0.02);
  EXPECT_NEAR(rhos[2][1], 0.000, 0.02);

  auto eq = evol.EquilibriumPopulationLowMut(rhos);
  // eq: [0.302589, 0.0434844, 0.653927]
  EXPECT_NEAR(eq[0], 0.30, 0.02);
  EXPECT_NEAR(eq[1], 0.04, 0.02);
  EXPECT_NEAR(eq[2], 0.66, 0.02);
}

TEST(EvolPrivRepGameAllCAllD, L1) {
  EvolPrivRepGame::SimulationParameters params;
  params.n_init = 1e5;
  params.n_steps = 1e5;

  EvolPrivRepGameAllCAllD evol(50, params, 5.0, 1.0);

  auto selfc_rho_eq = evol.EquilibriumCoopLevelAllCAllD(Norm::L1());
  double self_cooperation_level = std::get<0>(selfc_rho_eq);
  auto rhos = std::get<1>(selfc_rho_eq);
  auto eq = std::get<2>(selfc_rho_eq);

  // IC(self_cooperation_level, rhos, eq);
  // ic| self_cooperation_level: 0.90224
  // ic| rhos: [
  // [0, 0.0933583, 4.4365e-29],
  // [0.0120976, 0, 0.670243],
  // [0.0453502, 2.35795e-24, 0]
  // ]
  // ic| eq: [0.316563, 0.0433123, 0.640125]

  EXPECT_NEAR(self_cooperation_level, 0.90, 0.02);
  EXPECT_NEAR(rhos[0][1], 0.097, 0.02);
  EXPECT_NEAR(rhos[0][2], 0.000, 0.02);
  EXPECT_NEAR(rhos[1][0], 0.012, 0.02);
  EXPECT_NEAR(rhos[2][0], 0.043, 0.02);
  EXPECT_NEAR(eq[0], 0.30, 0.02);
  EXPECT_NEAR(eq[1], 0.04, 0.02);
  EXPECT_NEAR(eq[2], 0.66, 0.02);
}

TEST(EvolPrivRepGameFiniteMutationRateAllCAllD, L1) {
  EvolPrivRepGame::SimulationParameters params;
  params.n_init = 1e4;
  params.n_steps = 1e4;

  Norm L1 = Norm::L1();
  EvolPrivRepGameFiniteMutationRateAllCAllD evol(50, L1, params);
  double mu = 1.0e-2;
  double benefit = 5.0;
  auto result = evol.CalculateEquilibrium(benefit, 1.0, mu);
  EXPECT_NEAR(result.OverallCooperationLevel(), 0.74, 0.02);
  EXPECT_NEAR(result.OverallAbundances()[0], 0.47, 0.02);
  EXPECT_NEAR(result.OverallAbundances()[1], 0.32, 0.02);
  EXPECT_NEAR(result.OverallAbundances()[2], 0.21, 0.02);
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

std::vector<Norm> LocalMutants(const Norm& norm) {
  std::vector<Norm> local_mutants;

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
    local_mutants.push_back(mutant);
  }

  return local_mutants;
}

void PrintSelectionMutationEquilibriumAllCAllD(const Norm& norm, const SimulationParams& params, bool check_local_mutants = false) {
  auto start = std::chrono::high_resolution_clock::now();

  PrivateRepGame prg({{norm, params.N}}, params.seed);
  prg.Update(params.n_init, params.q, params.mu_percept, false);
  prg.ResetCounts();
  prg.Update(params.n_steps, params.q, params.mu_percept, true);
  IC( prg.NormAverageReputation(), prg.NormCooperationLevels());

  EvolPrivRepGame::SimulationParameters evo_params(params.n_init, params.n_steps, params.q, params.mu_percept, params.seed);
  EvolPrivRepGameAllCAllD evol(params.N, evo_params, params.benefit, params.beta);
  auto res = evol.EquilibriumCoopLevelAllCAllD(norm);
  auto rho = std::get<1>(res);
  auto eq = std::get<2>(res);
  IC(rho, eq);

  if (check_local_mutants) {
    CompareWithLocalMutants(norm, params);
  }

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cerr << "Elapsed time: " << elapsed.count() << " s\n";
}

void PrintESSCheck(const Norm& norm, const SimulationParams& params, bool check_local_mutants = false) {
  auto start = std::chrono::high_resolution_clock::now();

  auto compare_payoffs = [&params](const Norm& resident, const Norm& mutant) -> std::pair<double,double> {
    PrivateRepGame prg({{resident, params.N-1}, {mutant, 1}}, 123456789);
    prg.Update(params.n_init, params.q, params.mu_percept, false);
    prg.ResetCounts();
    prg.Update(params.n_steps, params.q, params.mu_percept, true);
    auto c_levels = prg.NormCooperationLevels();

    // pi_res = b[ (N-2)p_{res->res} + p_{mut->res} ] - c[ (N-2)p_{res->res} + p_{res->mut} ]
    // pi_mut = b[ (N-1)p_{res->mut}                ] - c[ (N-1)p_{mut->res}                ]
    //
    // delta_pi = pi_res - pi_mut = b[ (N-2)p_{res->res} + p_{mut->res} - (N-1)p_{res->mut} ] - c[ (N-2)p_{res->res} + p_{res->mut} - (N-1)p_{mut->res} ]
    //          = b P1 - c P2
    // if P1 >= 0; b/c > P2/P1
    // if P1 < 0; b/c < P2/P1
    double p_res_res = c_levels[0][0];
    double p_res_mut = c_levels[0][1];
    double p_mut_res = c_levels[1][0];
    double P1 = (params.N-2)*p_res_res + p_mut_res - (params.N-1)*p_res_mut;
    double P2 = (params.N-2)*p_res_res + p_res_mut - (params.N-1)*p_mut_res;
    if (P1 >= 0.0) {
      return {P2/P1, std::numeric_limits<double>::infinity()};
    }
    else {
      return {1.0, P2/P1};
    }
  };

  std::vector<Norm> mutants = {Norm::AllC(), Norm::AllD()};
  if (check_local_mutants) {
    auto locals = LocalMutants(norm);
    mutants.insert(mutants.end(), locals.begin(), locals.end());
  }

  std::cerr << "Calculating ESS b/c ranges:" << std::endl;
  double bc_min = 1.0, bc_max = std::numeric_limits<double>::infinity();
  for (const auto& mutant : mutants) {
    auto payoffs = compare_payoffs(norm, mutant);
    std::cerr << "0x" << std::setfill('0') << std::setw(5) << std::hex << mutant.ID() << std::resetiosflags(std::ios_base::fmtflags(-1))
              << " " << payoffs.first << ' ' << payoffs.second << std::endl;
    if (payoffs.first > bc_min) {
      bc_min = payoffs.first;
    }
    if (payoffs.second < bc_max) {
      bc_max = payoffs.second;
    }
  }
  std::cerr << "ESS b/c range: " << bc_min << " : " << bc_max << std::endl;

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

  std::vector<Norm> mutants = {Norm::AllC(), Norm::AllD()};
  auto local_mutants = LocalMutants(norm);
  mutants.insert(mutants.end(), local_mutants.begin(), local_mutants.end());
  double min_eq = 1.0;
  Norm min_eq_norm = Norm::AllC();
  for (const auto& mutant : mutants) {
    // std::cout << mutant.Inspect();
    EvolPrivRepGame evol(params.N, {norm, mutant}, evo_params);
    auto rhos = evol.FixationProbabilities(params.benefit, params.beta);
    auto eq = evol.EquilibriumPopulationLowMut(rhos);
    if (eq[0] < min_eq) {
      min_eq = eq[0];
      min_eq_norm = mutant;
    }
    std::cerr << "0x" << std::setfill('0') << std::setw(5) << std::hex << mutant.ID() << std::resetiosflags(std::ios_base::fmtflags(-1))
              << " " << eq[0] << std::endl;
  }
  std::cerr << "Most risky mutant:" << min_eq << std::endl;
  std::cerr << norm.InspectComparison(min_eq_norm);
}

/*
int main(int argc, char *argv[]) {

  std::vector<std::string> args;
  bool check_local_mutants = false;
  bool swap_gb = false;
  nlohmann::json j = nlohmann::json::object();
  // -j param.json : set parameters used for evolutionary simulation by json file
  // -l : check local mutants
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
    else if (std::string(argv[i]) == "-s") {
      swap_gb = true;
    }
    else {
      args.emplace_back(argv[i]);
    }
  }

  if (args.empty()) {
    test_RandomNorm();
    test_LeadingEight();
    test_SelectionMutationEquilibrium();
    test_SelectionMutationEquilibrium2();
    // test_SelectionMutationEquilibriumFiniteMu();
  }
  else if (args.size() == 1) {
    Norm n = Norm::ParseNormString(args.at(0), swap_gb);
    std::cout << n.Inspect();
    SimulationParams params = j.get<SimulationParams>();
    std::cout << nlohmann::json(params).dump(2) << std::endl;
    PrintSelectionMutationEquilibriumAllCAllD(n, params, check_local_mutants);
    PrintESSCheck(n, params, check_local_mutants);
  }
  else if (args.size() == 2) {  // if two arguments are given, direct competition between two norms are shown
    Norm n1 = Norm::ParseNormString(args.at(0), swap_gb);
    Norm n2 = Norm::ParseNormString(args.at(1), swap_gb);
    SimulationParams params = j.get<SimulationParams>();
    std::cout << nlohmann::json(params).dump(2) << std::endl;
    PrintCompetition(n1, n2, params);
  }

  return 0;
}
 */