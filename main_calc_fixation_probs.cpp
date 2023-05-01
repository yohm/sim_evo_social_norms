#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <tuple>
#include <chrono>
#include <mpi.h>
#include <nlohmann/json.hpp>
#include "Norm.hpp"
#include "PrivRepGame.hpp"


constexpr Reputation B = Reputation::B, G = Reputation::G;
constexpr Action C = Action::C, D = Action::D;

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

template <typename T>
class vector2d {
public:
  vector2d(size_t n_rows, size_t n_cols, T init) : _data(n_rows*n_cols, init), n_rows(n_rows), n_cols(n_cols) {};
  T& operator()(size_t i, size_t j) { return _data[i*n_cols+j]; }
  const T& operator()(size_t i, size_t j) const { return _data[i*n_cols+j]; }
  size_t Rows() const { return n_rows; }
  size_t Cols() const { return n_cols; }
  size_t size() const { return _data.size(); }
  T* data() { return _data.data(); }
  std::vector<T> _data;
  size_t n_rows;
  size_t n_cols;
};

double SelfCoopLevel(const Norm& norm, const SimulationParams& params) {
  PrivateRepGame prg({{norm, params.N}}, params.seed);
  prg.Update(params.n_init, params.q, params.mu_percept, false);
  prg.ResetCounts();
  prg.Update(params.n_steps, params.q, params.mu_percept, true);
  return prg.NormCooperationLevels()[0][0];
}

std::tuple<std::vector<int>, std::vector<double>, vector2d<double>> CalculateFixationProbsThirdOrderWithoutR2(const SimulationParams& params) {

  // loop over Norm
  constexpr size_t N_NORMS = 4096;
  std::vector<int> norm_ids(N_NORMS, 0);
  std::vector<double> self_coop_levels(N_NORMS, 0.0);
  vector2d<double> p_fix(N_NORMS, N_NORMS, 0.0);

  auto idx = [] (const Norm& n) { return n.IDwithoutR2(); };

  std::vector<Norm> unique_norms;
  std::vector<size_t> norm_index(N_NORMS, 0);
  for (int i = 0; i < N_NORMS; i++) {
    Norm norm = Norm::ConstructFromIDwithoutR2(i);
    norm_ids[i] = norm.ID();

    if ( i < idx(norm.SwapGB()) ) {
      norm_index[i] = idx(norm.SwapGB());
    }
    else {
      norm_index[i] = i;
      unique_norms.push_back(norm);
    }
  }

  EvolPrivRepGame::SimulationParameters evoparams({params.n_init, params.n_steps, params.q, params.mu_percept, params.seed});

  int my_rank = 0, num_procs = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  for (size_t i = 0; i < unique_norms.size(); i++) {
    if (i % num_procs != my_rank) continue;
    const Norm& n1 = unique_norms[i];
    double pc = SelfCoopLevel(n1, params);
    self_coop_levels[idx(n1)] = pc;
    p_fix(idx(n1),idx(n1)) = 1.0 / params.N;
  }

  std::vector<std::array<size_t,2>> ij_pairs{};
  for (size_t i = 0; i < unique_norms.size(); i++) {
    for (size_t j = i+1; j < unique_norms.size(); j++) {
      ij_pairs.push_back({i, j});
    }
  }

  // loop over ij_pairs
  for (size_t t=0; t < ij_pairs.size(); t++) {
    if (t % num_procs != my_rank) continue;
    size_t i = ij_pairs[t][0];
    size_t j = ij_pairs[t][1];
    if (t % 10'000 == 0) {
      std::cerr << "t / t_max: " << t << " / " << ij_pairs.size() << std::endl;
    }
    //for (const auto& [i,j] : ij_pairs) {
    const Norm& n1 = unique_norms[i];
    const Norm& n2 = unique_norms[j];
    EvolPrivRepGame evol(params.N, std::vector<Norm>({n1, n2}), evoparams);
    auto rhos = evol.FixationProbabilities(params.benefit, params.beta);
    p_fix(idx(n1),idx(n2)) = rhos[0][1];
    p_fix(idx(n2),idx(n1)) = rhos[1][0];
  }

  // take the sum of p_fix using MPI
  MPI_Allreduce(MPI_IN_PLACE, self_coop_levels.data(), self_coop_levels.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, p_fix.data(), p_fix.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  // calculate non-unique-norms
  for (size_t i = 0; i < N_NORMS; i++) {
    size_t ni = norm_index[i];
    self_coop_levels[i] = self_coop_levels[ni];
    for (size_t j = 0; j < N_NORMS; j++) {
      size_t nj = norm_index[j];
      p_fix(i,j) = p_fix(ni,nj);
    }
  }

  return std::make_tuple(norm_ids, self_coop_levels, p_fix);
}

std::tuple<std::vector<int>, std::vector<double>, vector2d<double>> CalculateFixationProbs(const SimulationParams& params, const std::vector<Norm>& norms) {

  auto idx = [&norms] (const Norm& n)->int {
    return std::distance(norms.begin(), std::find(norms.begin(), norms.end(), n));
  };

  const size_t N_NORMS = norms.size();
  std::vector<int> norm_ids(N_NORMS, 0);
  std::vector<double> self_coop_levels(N_NORMS, 0.0);
  vector2d<double> p_fix(N_NORMS, N_NORMS, 0.0);

  std::vector<Norm> unique_norms;
  std::vector<size_t> norm_index(N_NORMS, 0);
  for (int i = 0; i < N_NORMS; i++) {
    Norm norm = norms[i];
    norm_ids[i] = norm.ID();

    Norm swapped = norm.SwapGB();
    int swap_idx = idx(swapped);
    if ( i < swap_idx && swap_idx < N_NORMS ) {
      norm_index[i] = swap_idx;
    }
    else {
      norm_index[i] = i;
      unique_norms.push_back(norm);
    }
  }

  EvolPrivRepGame::SimulationParameters evoparams({params.n_init, params.n_steps, params.q, params.mu_percept, params.seed});

  int my_rank = 0, num_procs = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  for (size_t i = 0; i < unique_norms.size(); i++) {
    if (i % num_procs != my_rank) continue;
    const Norm& n1 = unique_norms[i];
    double pc = SelfCoopLevel(n1, params);
    self_coop_levels[idx(n1)] = pc;
    p_fix(idx(n1),idx(n1)) = 1.0 / params.N;
  }

  std::vector<std::array<size_t,2>> ij_pairs{};
  for (size_t i = 0; i < unique_norms.size(); i++) {
    for (size_t j = i+1; j < unique_norms.size(); j++) {
      ij_pairs.push_back({i, j});
    }
  }

  // loop over ij_pairs
  for (size_t t=0; t < ij_pairs.size(); t++) {
    if (t % num_procs != my_rank) continue;
    size_t i = ij_pairs[t][0];
    size_t j = ij_pairs[t][1];
    if (t % 10'000 == 0) {
      std::cerr << "t / t_max: " << t << " / " << ij_pairs.size() << std::endl;
    }
    //for (const auto& [i,j] : ij_pairs) {
    const Norm& n1 = unique_norms[i];
    const Norm& n2 = unique_norms[j];
    EvolPrivRepGame evol(params.N, std::vector<Norm>({n1, n2}), evoparams);
    auto rhos = evol.FixationProbabilities(params.benefit, params.beta);
    p_fix(idx(n1),idx(n2)) = rhos[0][1];
    p_fix(idx(n2),idx(n1)) = rhos[1][0];
  }

  // take the sum of p_fix using MPI
  MPI_Allreduce(MPI_IN_PLACE, self_coop_levels.data(), self_coop_levels.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, p_fix.data(), p_fix.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  // calculate non-unique-norms
  for (size_t i = 0; i < N_NORMS; i++) {
    size_t ni = norm_index[i];
    self_coop_levels[i] = self_coop_levels[ni];
    for (size_t j = 0; j < N_NORMS; j++) {
      size_t nj = norm_index[j];
      p_fix(i,j) = p_fix(ni,nj);
    }
  }

  return std::make_tuple(norm_ids, self_coop_levels, p_fix);
}


void WriteInMsgpack(const std::string& filepath, const SimulationParams& params, const std::vector<int>& norm_ids, const std::vector<double>& self_coop_levels, const vector2d<double>& p_fix) {
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

void PrintFixationProbsInText(std::ofstream& out, const std::vector<int>& norm_ids, const std::vector<double>& self_coop_levels, const vector2d<double>& p_fix) {
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
  SimulationParams params = j.get<SimulationParams>();
  if (my_rank == 0) std::cerr << "params: " << nlohmann::json(params) << std::endl;

  // measure elapsed time
  auto start = std::chrono::system_clock::now();

  // second-order norms
  // std::vector<Norm> norms;
  // for (size_t i = 0; i < 4096; i++) {
  //   Norm norm = Norm::ConstructFromIDwithoutR2(i);
  //   if (norm.IsSecondOrder()) {
  //     norms.push_back(norm);
  //   }
  // }

  // three-species systems
  std::vector<Norm> norms = {Norm::L1(), Norm::AllC(), Norm::AllD()};

  const auto result = CalculateFixationProbs(params, norms);
  // auto result = CalculateFixationProbsThirdOrderWithoutR2(params);

  std::vector<int> norm_ids = std::get<0>(result);
  std::vector<double> self_coop_levels = std::get<1>(result);
  vector2d<double> p_fix = std::get<2>(result);

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

