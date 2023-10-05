#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <chrono>
#include <random>
#include <nlohmann/json.hpp>
#include <icecream.hpp>
#include "Norm.hpp"


class GroupedEvo {
  public:
  GroupedEvo(const std::vector<size_t> &norm_ids, const std::vector<std::vector<double>> &p_fix,
             const std::vector<double> &self_coop_levels) :
    norm_ids(norm_ids), p_fix(p_fix), self_coop_levels(self_coop_levels) {
    N_NORMS = norm_ids.size();
    assert(p_fix.size() == N_NORMS);
    assert(self_coop_levels.size() == N_NORMS);
  }

  size_t N_NORMS;
  std::vector<size_t> norm_ids;
  std::vector<std::vector<double>> p_fix;
  std::vector<double> self_coop_levels;

  static double InterGroupImitationProb(double pc_res, double pc_mut, double benefit, double sigma_out) {
    double pi_res = (benefit - 1.0) * pc_res;
    double pi_mut = (benefit - 1.0) * pc_mut;
    // f_{A\to B} = { 1 + \exp[ \sigma_out (s_A - s_B) ] }^{-1}
    return 1.0 / (1.0 + std::exp(sigma_out * (pi_res - pi_mut)));
  }

  double CalcAlpha(size_t i, size_t j, double benefit, double sigma_out) const {
    double p_plus = InterGroupImitationProb(self_coop_levels[i], self_coop_levels[j], benefit, sigma_out) * p_fix[i][j];
    double p_minus =
      InterGroupImitationProb(self_coop_levels[j], self_coop_levels[i], benefit, sigma_out) * p_fix[j][i];
    return p_plus - p_minus;
  }
};


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
    std::cout << "incoming flow from " << ni << " : " << i << std::endl;
    for (size_t j = 0; j < 10; ++j) {
      auto norm = Norm::ConstructFromID(alpha_ji[j].first);
      std::string name = norm.GetName();
      if (name.empty()) { name = "other";}
      std::cout << alpha_ji[j].first << ' ' << alpha_ji[j].second << ' ' << name << std::endl;
    }

    // print out the 10 smallest alpha_ji
    std::cout << "outgoing flow to " << ni << " : " << i << std::endl;
    for (size_t j = N_NORMS - 10; j < N_NORMS; ++j) {
      auto norm = Norm::ConstructFromID(alpha_ji[j].first);
      std::string name = norm.GetName();
      if (name.empty()) { name = "other";}
      std::cout << alpha_ji[j].first << ' ' << alpha_ji[j].second << ' ' << name << std::endl;
    }
  }
  else if (argc == 4) {
    size_t ni = std::stoul(argv[2]);
    size_t nj = std::stoul(argv[3]);
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

  }

  return 0;
}
