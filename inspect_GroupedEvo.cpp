#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <chrono>
#include <random>
#include <nlohmann/json.hpp>
#include <icecream.hpp>
#include "Norm.hpp"
#include "GroupedEvo.hpp"


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
  // calculate stationary state by ODE for evolution in group-structured population

  using json = nlohmann::json;

  if (argc != 3 && argc != 4) {
    std::cerr << "Usage: " << argv[0] << " <input_msgpack_file> <norm_id1> [norm_id2]" << std::endl;
    return 1;
  }

  // load fixation probabilities from the input file
  json j_in = LoadMsgpackFile(argv[1]);
  std::cerr << "json loaded" << std::endl;

  double sigma_out = j_in["sigma_in"].get<double>();
  double benefit = j_in["benefit"].get<double>();

  std::vector<size_t> norm_ids;
  for (auto& norm_id_j : j_in["norm_ids"]) {
    int id = norm_id_j.get<int>();
    norm_ids.emplace_back(id);
  }
  std::cerr << "  # of norms: " << norm_ids.size() << std::endl;

  const size_t N_NORMS = norm_ids.size();
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

  GroupedEvo ge(norm_ids, p_fix, self_coop_levels);

  if (argc == 3) {
    // size_t ni = std::stoul(argv[2]);
    Norm norm = Norm::ParseNormString(argv[2]);
    if (norm.IDwithoutR2() < norm.SwapGB().IDwithoutR2()) {
      norm = norm.SwapGB();
    }
    // print largest
    size_t ni = norm.ID();
    // find index from norm_ids
    size_t i = 0;
    for (; i < N_NORMS; ++i) {
      if (norm_ids[i] == ni) {
        break;
      }
    }
    if (i == N_NORMS) {
      std::cerr << "norm_id " << ni << " not found" << std::endl;
      return 1;
    }

    std::vector<std::pair<size_t, double>> alpha_ji;
    for (size_t j = 0; j < N_NORMS; ++j) {
      alpha_ji.emplace_back(norm_ids[j], ge.CalcAlpha(j, i, benefit, sigma_out) );
    }
    // sort alpha in descending order
    std::sort(alpha_ji.begin(), alpha_ji.end(), [](auto& a, auto& b) {
      return a.second > b.second;
    });

    // print out the 10 largest alpha_ji
    std::cout << "incoming flow (positive alpha_ji) from " << ni << " : " << i << std::endl;
    for (size_t j = 0; j < 10; ++j) {
      if (alpha_ji[j].second < 0.0) { break; }
      auto norm_j = Norm::ConstructFromID(alpha_ji[j].first);
      std::string name = norm_j.GetName();
      if (name.empty()) { name = "other";}
      std::cout << alpha_ji[j].first << ' ' << alpha_ji[j].second << ' ' << name << std::endl;
    }

    // print out the 10 smallest alpha_ji
    std::cout << "outgoing flow (negative alpha_ji) to " << ni << " : " << i << std::endl;
    for (size_t j = N_NORMS - 10; j < N_NORMS; ++j) {
      if (alpha_ji[j].second > 0.0) { break; }
      auto norm_j = Norm::ConstructFromID(alpha_ji[j].first);
      std::string name = norm_j.GetName();
      if (name.empty()) { name = "other";}
      std::cout << alpha_ji[j].first << ' ' << alpha_ji[j].second << ' ' << name << std::endl;
    }

    // mutational outflow  \sum_q p_fix[i][j]
    std::cout << "mutational outflow from " << ni << " : " << i << std::endl;
    std::cout << ge.MutationOutFlow(i) << std::endl;
  }
  else if (argc == 4) {
    Norm norm1 = Norm::ParseNormString(argv[2]);
    if (norm1.IDwithoutR2() < norm1.SwapGB().IDwithoutR2()) {
      norm1 = norm1.SwapGB();
    }
    Norm norm2 = Norm::ParseNormString(argv[3]);
    if (norm2.IDwithoutR2() < norm2.SwapGB().IDwithoutR2()) {
      norm2 = norm2.SwapGB();
    }
    size_t ni = norm1.ID();
    size_t nj = norm2.ID();
    // find index from norm_ids
    size_t i = 0;
    for (; i < N_NORMS; ++i) {
      if (norm_ids[i] == ni) {
        break;
      }
    }
    if (i == N_NORMS) {
      std::cerr << "norm_id " << ni << " not found" << std::endl;
      return 1;
    }
    size_t j = 0;
    for (; j < N_NORMS; ++j) {
      if (norm_ids[j] == nj) {
        break;
      }
    }
    if (j == N_NORMS) {
      std::cerr << "norm_id " << nj << " not found" << std::endl;
      return 1;
    }

    std::cout << "flow " << ni << " -> " << nj << " : " << ge.CalcAlpha(i, j, benefit, sigma_out) << std::endl;
    std::cout << "flow " << nj << " -> " << ni << " : " << ge.CalcAlpha(j, i, benefit, sigma_out) << std::endl;
    std::cout << "p_fix " << ni << " -> " << nj << " : " << p_fix[i][j] << std::endl;
    std::cout << "p_fix " << nj << " -> " << ni << " : " << p_fix[j][i] << std::endl;
  }

  return 0;
}
