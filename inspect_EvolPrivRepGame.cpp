#include <iostream>
#include <fstream>
#include <chrono>
#include <regex>
#include <icecream.hpp>
#include <nlohmann/json.hpp>
#include "EvolPrivRepGame.hpp"


void PrintSelectionMutationEquilibriumAllCAllD(const Norm& norm, const EvolPrivRepGame::Parameters& params, double benefit, double beta) {
  auto start = std::chrono::high_resolution_clock::now();

  std::cerr << "Norm: " << norm.GetName() << std::endl;
  std::cerr << "Parameteres:" << nlohmann::json(params).dump(2)  << std::endl;
  std::cerr << "benefit: " << benefit << std::endl;
  std::cerr << "beta: " << beta << std::endl;

  auto selfc_rho_eq = EvolPrivRepGame::EquilibriumCoopLevelAllCAllD(norm, params, benefit, beta);

  double self_cooperation_level = std::get<0>(selfc_rho_eq);
  auto rhos = std::get<1>(selfc_rho_eq);
  auto eq = std::get<2>(selfc_rho_eq);
  double eq_cooperation_level = self_cooperation_level * eq[0] + 1.0 * eq[1] + 0.0 * eq[2];
  IC(self_cooperation_level, rhos, eq, eq_cooperation_level);

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cerr << "Elapsed time: " << elapsed.count() << " s\n";
}


void PrintCompetition(const Norm& n1, const Norm& n2, const EvolPrivRepGame::Parameters& params, double benefit, double beta) {
  auto start = std::chrono::high_resolution_clock::now();

  int omp_threads = omp_get_max_threads();
  std::cerr << "Running with " << omp_threads << " threads" << std::endl;
  std::cerr << "Parameteres:" << nlohmann::json(params).dump(2)  << std::endl;
  std::cerr << "benefit: " << benefit << std::endl;
  std::cerr << "beta: " << beta << std::endl;

  auto rhoij_rhoji_pii_pij = EvolPrivRepGame::FixationProbabilityAndPayoff(n1, n2, params, benefit, beta);
  // auto rhos = EvolPrivRepGame::FixationProbabilities({n1, n2}, params, benefit, beta);

  std::vector<std::vector<double>> rhos = { {0.0, std::get<0>(rhoij_rhoji_pii_pij)}, {std::get<1>(rhoij_rhoji_pii_pij), 0.0} };

  auto eq = EvolPrivRepGame::EquilibriumPopulationLowMut(rhos);
  std::cerr << "Transition probabilities between " << n1.GetName() << " vs " << n2.GetName() << std::endl;
  std::cerr << "  ---> : " << rhos[0][1] << std::endl;
  std::cerr << "  <--- : " << rhos[1][0] << std::endl;
  std::cerr << "Equilibrium population: " << eq[0] << " , " << eq[1] << std::endl;

  double pc_s1 = EvolPrivRepGame::MonomorphicCooperationLevel(n1, params);
  double pc_s2 = EvolPrivRepGame::MonomorphicCooperationLevel(n2, params);
  std::cerr << "Monomorphic cooperation levels: " << pc_s1 << " , " << pc_s2 << std::endl;

  auto pi_i = std::get<2>(rhoij_rhoji_pii_pij);
  auto pi_j = std::get<3>(rhoij_rhoji_pii_pij);
  std::cerr << "# num mutants l, pi_i[l], pi_j[l]" << std::endl;
  for (size_t l = 1; l < params.N; l++) {
    std::cout << l << " " << pi_i[l] << " " << pi_j[l] << std::endl;
  }

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cerr << "Elapsed time: " << elapsed.count() << " s\n";
}

int main(int argc, char *argv[]) {

  std::vector<std::string> args;
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
    else {
      args.emplace_back(argv[i]);
    }
  }

  EvolPrivRepGame::Parameters params = j.get<EvolPrivRepGame::Parameters>();
  // if j has a key "benefit" and "beta, use these as a benefit value
  double benefit = 5.0, beta = 1.0;
  if (j.contains("benefit")) {
    benefit = j["benefit"].get<double>();
  }
  if (j.contains("beta")) {
    beta = j["beta"].get<double>();
  }

  if (args.size() == 1) {
    Norm s = Norm::ParseNormString(args.at(0));
    PrintSelectionMutationEquilibriumAllCAllD(s, params, benefit, beta);
  }
  else if (args.size() == 2) {  // if two arguments are given, direct competition between two norms are shown
    Norm s1 = Norm::ParseNormString(args.at(0));
    Norm s2 = Norm::ParseNormString(args.at(1));
    PrintCompetition(s1, s2, params, benefit, beta);
  }
  else {
    std::cerr << "Usage: " << argv[0] << " [-j param.json] norm_string" << std::endl;
    std::cerr << "       " << argv[0] << " [-j param.json] norm_string1 norm_string2" << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << "  -j param.json : set parameters used for evolutionary simulation by json file" << std::endl;
    std::cerr << "  norm_string : string representation of a norm" << std::endl;
    std::cerr << "Default parameters:" << std::endl;
    std::cerr << "  " << nlohmann::json(EvolPrivRepGame::Parameters{}).dump(2) << std::endl;
    return 1;
  }

  return 0;
}