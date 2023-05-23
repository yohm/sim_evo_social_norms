#include <iostream>
#include <fstream>
#include <chrono>
#include <regex>
#include <icecream.hpp>
#include <nlohmann/json.hpp>
#include "PrivRepGame.hpp"
#include "Parameters.hpp"


void CompareWithLocalMutants(const Norm& norm, const Parameters& params);

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

void PrintSelectionMutationEquilibriumAllCAllD(const Norm& norm, const Parameters& params, bool check_local_mutants = false) {
  auto start = std::chrono::high_resolution_clock::now();

  PrivateRepGame prg({{norm, params.N}}, params.seed);
  prg.Update(params.n_init, params.q, params.mu_percept, false);
  prg.ResetCounts();
  prg.Update(params.n_steps, params.q, params.mu_percept, true);
  IC( prg.NormAverageReputation(), prg.NormCooperationLevels());

  auto evo_params = params.ToEvolParams();
  EvolPrivRepGameAllCAllD evol(evo_params, params.benefit, params.beta);
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

void PrintESSCheck(const Norm& norm, const Parameters& params, bool check_local_mutants = false) {
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

void PrintCompetition(const Norm& n1, const Norm& n2, const Parameters& params) {
  auto start = std::chrono::high_resolution_clock::now();

  EvolPrivRepGame::SimulationParameters evo_params = params.ToEvolParams();

  EvolPrivRepGame evol(evo_params);
  auto fixs = evol.FixationProbabilities({n1, n2}, params.benefit, params.beta);
  auto eq = evol.EquilibriumPopulationLowMut(fixs);
  IC( fixs, eq );

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cerr << "Elapsed time: " << elapsed.count() << " s\n";
}

void CompareWithLocalMutants(const Norm& norm, const Parameters& params) {
  EvolPrivRepGame::SimulationParameters evo_params = params.ToEvolParams();

  std::vector<Norm> mutants = {Norm::AllC(), Norm::AllD()};
  auto local_mutants = LocalMutants(norm);
  mutants.insert(mutants.end(), local_mutants.begin(), local_mutants.end());
  double min_eq = 1.0;
  Norm min_eq_norm = Norm::AllC();
  for (const auto& mutant : mutants) {
    // std::cout << mutant.Inspect();
    EvolPrivRepGame evol(evo_params);
    auto rhos = evol.FixationProbabilities({norm, mutant}, params.benefit, params.beta);
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

  if (args.size() == 1) {
    Norm n = Norm::ParseNormString(args.at(0), swap_gb);
    std::cout << n.Inspect();
    Parameters params = j.get<Parameters>();
    std::cout << nlohmann::json(params).dump(2) << std::endl;
    PrintSelectionMutationEquilibriumAllCAllD(n, params, check_local_mutants);
    PrintESSCheck(n, params, check_local_mutants);
  }
  else if (args.size() == 2) {  // if two arguments are given, direct competition between two norms are shown
    Norm n1 = Norm::ParseNormString(args.at(0), swap_gb);
    Norm n2 = Norm::ParseNormString(args.at(1), swap_gb);
    Parameters params = j.get<Parameters>();
    std::cout << nlohmann::json(params).dump(2) << std::endl;
    PrintCompetition(n1, n2, params);
  }
  else {
    std::cerr << "Usage: " << argv[0] << " [-j param.json] [-l] [-s] norm_string" << std::endl;
    std::cerr << "       " << argv[0] << " [-j param.json] [-l] [-s] norm_string1 norm_string2" << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << "  -j param.json : set parameters used for evolutionary simulation by json file" << std::endl;
    std::cerr << "  -l : check local mutants" << std::endl;
    std::cerr << "  -s : swap good and bad" << std::endl;
    std::cerr << "  norm_string : string representation of a norm" << std::endl;
    std::cerr << "Default parameters:" << std::endl;
    std::cerr << "  " << nlohmann::json(Parameters()).dump(2) << std::endl;
    return 1;
  }

  return 0;
}