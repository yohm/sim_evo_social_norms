#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <chrono>
#include <nlohmann/json.hpp>
#include "Vector2d.hpp"
#include "Norm.hpp"
#include "PrivRepGame.hpp"
#include "Parameters.hpp"


void PrintFixationProbParamDep(const Parameters& params, const Norm& resident, const Norm& mutant, std::ostream& out) {
  std::vector<double> sigma_array = {0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0};
  std::vector<std::pair<double,double>> benefit_sigma_vec;
  for (double sigma : sigma_array) {
    for (int b_int = 105; ; b_int += 5) {
      double benefit = static_cast<double>(b_int) / 100.0;
      constexpr double benefit_max = 6.0;
      if (benefit > benefit_max) break;
      benefit_sigma_vec.emplace_back(benefit, sigma);
    }
  }

  auto evoparams = params.ToEvolParams();
  EvolPrivRepGame evol(evoparams);
  std::vector<std::pair<double,double>> fixs = evol.FixationProbabilityBatch(resident, mutant, benefit_sigma_vec);

  for (size_t i = 0; i < benefit_sigma_vec.size(); i++) {
    out << benefit_sigma_vec[i].first << ' ' << benefit_sigma_vec[i].second << ' ' << fixs[i].first << std::endl;
  }
}


int main(int argc, char* argv[]) {

  using namespace nlohmann;

  json j = json::object();
  std::vector<std::string> args;

  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == "-j" && i + 1 < argc) {
      std::ifstream fin(argv[++i]);
      // check if file exists
      if (fin) {
        fin >> j;
        fin.close();
      } else {
        std::istringstream iss(argv[i]);
        iss >> j;
      }
    }
    else {
      args.emplace_back(argv[i]);
    }
  }
  Parameters params = j.get<Parameters>();

  std::cerr << "params: " << nlohmann::json(params) << std::endl;

  if (args.size() != 2) {
    std::cerr << "Error! : two arguments are necessary" << std::endl;
    std::cerr << "Usage: " << args[0] << " <resident> <mutant> [-j param.json]" << std::endl;
    for (auto& arg : args) {
      std::cerr << "arg: " << arg << std::endl;
    }
    return 1;
  }

  std::cerr << "resident: " << args[0] << std::endl;
  std::cerr << "mutant: " << args[1] << std::endl;

  // measure elapsed time
  auto start = std::chrono::system_clock::now();

  auto resident = Norm::ParseNormString(args[0]);
  auto mutant = Norm::ParseNormString(args[1]);

  PrintFixationProbParamDep(params, resident, mutant, std::cout);

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cerr << "elapsed time: " << elapsed_seconds.count() << "s\n";

  return 0;
}

