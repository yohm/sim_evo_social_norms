#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <chrono>
#include <mpi.h>
#include <nlohmann/json.hpp>
#include "Vector2d.hpp"
#include "Norm.hpp"
#include "PrivRepGame.hpp"
#include "Parameters.hpp"


double SelfCoopLevel(const Norm& norm, const EvolPrivRepGame::SimulationParameters& evoparams) {
  EvolPrivRepGame evol(evoparams);
  return evol.SelfCooperationLevel(norm);
}

using p_fix_vec_t = std::vector< Vector2d<double> >;
void CalculateFixationProbs(const ParametersBatch& params, const std::vector<Norm>& norms, std::vector<double>& self_coop_levels, p_fix_vec_t& p_fix_vec) {

  const size_t NN = norms.size();   // number of norms
  self_coop_levels.clear();
  self_coop_levels.assign(NN, 0.0);
  size_t NP = params.benefit_sigma_in_vec().size();
  p_fix_vec.clear();
  p_fix_vec.assign(NP, Vector2d<double>(NN, NN, 0.0));

  EvolPrivRepGame::SimulationParameters evoparams = params.ToEvolParams();

  int my_rank = 0, num_procs = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  // calculate self-cooperation levels
  for (size_t i = 0; i < NN; i++) {
    if (i % num_procs != my_rank) continue;
    const Norm& n1 = norms[i];
    evoparams.seed += i;
    double pc = SelfCoopLevel(n1, evoparams);
    self_coop_levels[i] = pc;
    for (auto& p_fix : p_fix_vec) {
      p_fix(i, i) = 1.0 / params.N;
    }
  }

  std::vector<std::array<size_t,2>> ij_pairs{};
  for (size_t i = 0; i < NN; i++) {
    for (size_t j = i+1; j < NN; j++) {
      ij_pairs.push_back({i, j});
    }
  }

  // loop over ij_pairs
  for (size_t ij=0; ij < ij_pairs.size(); ij++) {
    if (ij % num_procs != my_rank) continue;
    size_t i = ij_pairs[ij][0];
    size_t j = ij_pairs[ij][1];
    if (ij % 10'000 == 0) {
      std::cerr << "ij / ij_max: " << ij << " / " << ij_pairs.size() << std::endl;
    }
    const Norm& n1 = norms[i];
    const Norm& n2 = norms[j];
    evoparams.seed += NN + ij * params.N;
    EvolPrivRepGame evol(evoparams);
    auto fs_vec = evol.FixationProbabilityBatch(n1, n2, params.benefit_sigma_in_vec());
    for (size_t n = 0; n < NP; n++) {
      p_fix_vec[n](i,j) = fs_vec[n].first;
      p_fix_vec[n](j,i) = fs_vec[n].second;
    }
  }

  // take the sum of p_fix using MPI
  MPI_Allreduce(MPI_IN_PLACE, self_coop_levels.data(), self_coop_levels.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for (size_t n = 0; n < NP; n++) {
    MPI_Allreduce(MPI_IN_PLACE, p_fix_vec[n].data(), p_fix_vec[n].size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
}


void WriteInMsgpack(const std::string& filepath, const Parameters& params, const std::vector<int>& norm_ids, const std::vector<double>& self_coop_levels, const Vector2d<double>& p_fix) {
  // convert to json
  nlohmann::json j_out = nlohmann::json::object();
  j_out["params"] = params;
  j_out["norm_ids"] = norm_ids;
  j_out["self_coop_levels"] = self_coop_levels;
  j_out["p_fix"] = p_fix._data;
  // write a messagepack into a binary file using json-library
  std::ofstream ofs(filepath, std::ios::binary);
  std::vector<std::uint8_t> v = nlohmann::json::to_msgpack(j_out);
  ofs.write(reinterpret_cast<const char*>(v.data()), v.size());
  ofs.close();
}

void PrintFixationProbsInText(std::ofstream& out, const std::vector<int>& norm_ids, const std::vector<double>& self_coop_levels, const Vector2d<double>& p_fix) {
  for (size_t i = 0; i < p_fix.Rows(); i++) {
    out << norm_ids[i] << " " << self_coop_levels[i] << " ";
    for (size_t j = 0; j < p_fix.Cols(); j++) {
      out << p_fix(i,j) << " ";
    }
    out << std::endl;
  }
}

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);

  int my_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  using namespace nlohmann;

  json j = json::object();
  bool debug_mode = false;
  bool leading_eight_only = false;
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
    } else if (std::string(argv[i]) == "-d") {
      debug_mode = true;
    } else if (std::string(argv[i]) == "-l8") {
      leading_eight_only = true;
    } else {
      std::cerr << "unknown option: " << argv[i] << std::endl;
      return 1;
    }
  }
  ParametersBatch params = j.get<ParametersBatch>();

  if (debug_mode) {
    // debug mode
    if (my_rank == 0) {
      std::cerr << "debug mode" << std::endl;
      std::cerr << "params: " << nlohmann::json(params) << std::endl;
    }
    MPI_Finalize();
    return 0;
  }

  if (my_rank == 0) std::cerr << "params: " << nlohmann::json(params) << std::endl;


  // measure elapsed time
  auto start = std::chrono::system_clock::now();

  std::vector<Norm> norms;
  if (leading_eight_only) {
    norms = {Norm::L1(), Norm::L2(), Norm::L3(), Norm::L4(), Norm::L5(), Norm::L6(), Norm::L7(), Norm::L8(), Norm::AllC(), Norm::AllD()};
  }
  else {
    norms = Norm::Deterministic3rdOrderWithoutR2Norms();
  }

  std::vector<double> self_coop_levels;
  p_fix_vec_t p_fix_vec;
  CalculateFixationProbs(params, norms, self_coop_levels, p_fix_vec);

  std::vector<int> norm_ids;
  for (const Norm& n: norms) {
    norm_ids.emplace_back(n.ID());
  }

  if (my_rank == 0) {
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cerr << "elapsed time: " << elapsed_seconds.count() << "s\n";

    for (size_t n = 0; n < p_fix_vec.size(); n++) {
      std::string filename = "fixation_probs_" + std::to_string(n) + ".msgpack";
      WriteInMsgpack(filename, params.ParameterAt(n), norm_ids, self_coop_levels, p_fix_vec[n]);
      // write txt file as well
      std::ofstream fout("fixation_probs_" + std::to_string(n) + ".dat");
      PrintFixationProbsInText(fout, norm_ids, self_coop_levels, p_fix_vec[n]);
    }
  }

  MPI_Finalize();
  return 0;
}

