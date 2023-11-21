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

  if (argc != 6) {
    std::cerr << "Usage: " << argv[0] << " <input_msgpack_file> <sigma_out> <relative mutation rate> <tmax> <dt>" << std::endl;
    return 1;
  }

  // load fixation probabilities from the input file
  json j_in = LoadMsgpackFile(argv[1]);
  std::cerr << "json loaded" << std::endl;

  double sigma_out = std::stod(argv[2]);
  double relative_mutation_rate = std::stod(argv[3]);
  size_t T_max = std::stoul(argv[4]);
  double dt = std::stod(argv[5]);

  std::vector<Norm> norms;
  std::vector<size_t> norm_ids;
  for (auto& norm_id_j : j_in["norm_ids"]) {
    int id = norm_id_j.get<int>();
    norm_ids.emplace_back(id);
    Norm n = Norm::ConstructFromID(id);
    norms.emplace_back(n);
  }
  std::cerr << "  # of norms: " << norms.size() << std::endl;

  std::ofstream norm_out("norms.txt");
  for (auto& n : norms) {
    std::string name = n.GetName();
    if (name.empty()) { name = '-'; }
    norm_out << n.ID() << ' ' << name << std::endl;
  }
  norm_out.close();

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
  double benefit = j_in["benefit"].get<double>();

  auto start = std::chrono::system_clock::now();

  GroupedEvo g_evo{norm_ids, p_fix, self_coop_levels};
  std::ofstream ts_out("timeseries.dat");
  auto stationary = g_evo.TimeEvolutionODE(benefit, sigma_out, relative_mutation_rate, T_max, dt, ts_out);

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cerr << "elapsed time: " << elapsed_seconds.count() << "s\n";


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
  std::ofstream fout("histo_norms.dat");
  for (auto& n : norms_to_measure) {
    fout << std::get<0>(n) << ' ' << std::get<1>(n) << ' ' << std::get<2>(n) << std::endl;
  }

  return 0;
}
