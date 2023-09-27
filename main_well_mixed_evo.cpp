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


int main(int argc, char* argv[]) {
  // run evolutionary simulation in homogeneous population

  using json = nlohmann::json;

  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <input_msgpack_file>" << std::endl;
    return 1;
  }

  // load fixation probabilities from the input file
  json j_in = LoadMsgpackFile(argv[1]);
  std::cerr << "json loaded" << std::endl;

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

  auto stationary = EvolPrivRepGame::EquilibriumPopulationLowMutPowerMethod(p_fix);

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
  for (auto& n : norms_to_measure) {
    std::cout << std::get<0>(n) << ' ' << std::get<1>(n) << ' ' << std::get<2>(n) << std::endl;
  }

  return 0;
}
