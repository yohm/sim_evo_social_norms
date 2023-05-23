#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <tuple>
#include <map>
#include <chrono>
#include <mpi.h>
#include <nlohmann/json.hpp>
#include "Vector2d.hpp"
#include "Norm.hpp"
#include "PrivRepGame.hpp"
#include "Parameters.hpp"


constexpr Reputation B = Reputation::B, G = Reputation::G;
constexpr Action C = Action::C, D = Action::D;

double SelfCoopLevel(const Norm& norm, const Parameters& params) {
  PrivateRepGame prg({{norm, params.N}}, params.seed);
  prg.Update(params.n_init, params.q, params.mu_percept, false);
  prg.ResetCounts();
  prg.Update(params.n_steps, params.q, params.mu_percept, true);
  return prg.NormCooperationLevels()[0][0];
}

std::pair<std::vector<double>, Vector2d<double>> CalculateFixationProbs(const Parameters& params, const std::vector<Norm>& norms) {

  const size_t NN = norms.size();   // number of norms
  std::vector<double> self_coop_levels(NN, 0.0);
  Vector2d<double> p_fix(NN, NN, 0.0);

  EvolPrivRepGame::SimulationParameters evoparams = params.ToEvolParams();

  int my_rank = 0, num_procs = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  // calculate self-cooperation levels
  for (size_t i = 0; i < NN; i++) {
    if (i % num_procs != my_rank) continue;
    const Norm& n1 = norms[i];
    double pc = SelfCoopLevel(n1, params);
    self_coop_levels[i] = pc;
    p_fix(i, i) = 1.0 / params.N;
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
    EvolPrivRepGame evol(evoparams);
    auto fs = evol.FixationProbability(n1, n2, params.benefit, params.beta);
    p_fix(i,j) = fs.first;
    p_fix(j,i) = fs.second;
  }

  // take the sum of p_fix using MPI
  MPI_Allreduce(MPI_IN_PLACE, self_coop_levels.data(), self_coop_levels.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, p_fix.data(), p_fix.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return std::make_pair(self_coop_levels, p_fix);
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
  // run evolutionary simulation in group-structured population
  // strategy space: deterministic strategies without R2 (~ 2000 strategies)

  MPI_Init(&argc, &argv);

  int my_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  using namespace nlohmann;

  std::vector<std::string> args;
  json j = json::object();
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
      std::cerr << "unknown option: " << argv[i] << std::endl;
      return 1;
    }
  }
  Parameters params = j.get<Parameters>();
  if (my_rank == 0) std::cerr << "params: " << nlohmann::json(params) << std::endl;

  // measure elapsed time
  auto start = std::chrono::system_clock::now();

  // std::vector<Norm> norms = {Norm::L1(), Norm::AllC(), Norm::AllD()};
  // std::vector<Norm> norms = Norm::Deterministic2ndOrderWithoutR2Norms();
  std::vector<Norm> norms = Norm::Deterministic3rdOrderWithoutR2Norms();
  const auto result = CalculateFixationProbs(params, norms);

  std::vector<int> norm_ids;
  for (const Norm& n: norms) {
    norm_ids.emplace_back(n.ID());
  }
  std::vector<double> self_coop_levels = result.first;
  Vector2d<double> p_fix = result.second;

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  if (my_rank == 0) {
    std::cerr << "elapsed time: " << elapsed_seconds.count() << "s\n";

    WriteInMsgpack("fixation_probs.msgpack", params, norm_ids, self_coop_levels, p_fix);
    std::ofstream fout("fixation_probs.dat");
    PrintFixationProbsInText(fout, norm_ids, self_coop_levels, p_fix);
  }

  MPI_Finalize();
  return 0;
}

