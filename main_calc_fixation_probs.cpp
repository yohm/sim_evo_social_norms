#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <chrono>
#include <mpi.h>
#include <nlohmann/json.hpp>
#include "Vector2d.hpp"
#include "Norm.hpp"
#include "EvolPrivRepGame.hpp"


double SelfCoopLevel(const Norm& norm, const EvolPrivRepGame::Parameters& evoparam) {
  return EvolPrivRepGame::MonomorphicCooperationLevel(norm, evoparam);
}

using p_fix_vec_t = std::vector< Vector2d<double> >;
using b_sigma_vec_t = std::vector< std::pair<double,double> >;
void CalculateFixationProbs(const EvolPrivRepGame::Parameters& evoparam, const b_sigma_vec_t& b_sigma_vec, const std::vector<Norm>& norms, std::vector<double>& self_coop_levels, p_fix_vec_t& p_fix_vec) {

  const size_t NN = norms.size();   // number of norms
  self_coop_levels.clear();
  self_coop_levels.assign(NN, 0.0);
  size_t NP = b_sigma_vec.size();
  p_fix_vec.clear();
  p_fix_vec.assign(NP, Vector2d<double>(NN, NN, 0.0));


  int my_rank = 0, num_procs = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  // calculate self-cooperation levels
  for (size_t i = 0; i < NN; i++) {
    if (i % num_procs != my_rank) continue;
    const Norm& n1 = norms[i];
    EvolPrivRepGame::Parameters params = evoparam;
    params.seed += i;
    double pc = SelfCoopLevel(n1, evoparam);
    self_coop_levels[i] = pc;
    for (auto& p_fix : p_fix_vec) {
      p_fix(i, i) = 1.0 / static_cast<double>(params.N);
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
    EvolPrivRepGame::Parameters params = evoparam;
    params.seed += NN + ij * params.N;
    auto fs_vec = EvolPrivRepGame::FixationProbabilityBatch(n1, n2, params, b_sigma_vec);
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


void WriteInMsgpack(const std::string& filepath, const EvolPrivRepGame::Parameters& evoparam, double benefit, double sigma_in, const std::string& norm_set_name, const std::vector<int>& norm_ids, const std::vector<double>& self_coop_levels, const Vector2d<double>& p_fix) {
  // convert to json
  nlohmann::json j_out = nlohmann::json::object();
  j_out["evoparam"] = evoparam;
  j_out["benefit"] = benefit;
  j_out["sigma_in"] = sigma_in;
  j_out["norm_set"] = norm_set_name;
  j_out["norm_ids"] = norm_ids;
  j_out["self_coop_levels"] = self_coop_levels;
  j_out["p_fix"] = p_fix._data;
  // write a messagepack into a binary file using json-library
  std::ofstream ofs(filepath, std::ios::binary);
  std::vector<std::uint8_t> v = nlohmann::json::to_msgpack(j_out);
  ofs.write(reinterpret_cast<const char*>(v.data()), v.size());
  ofs.close();
}

void PrintFixationProbsInText(std::ofstream& out, const std::string& norm_set_name, const std::vector<int>& norm_ids, const std::vector<double>& self_coop_levels, const Vector2d<double>& p_fix) {
  out << "# norm_set: " << norm_set_name << std::endl;
  for (size_t i = 0; i < p_fix.Rows(); i++) {
    out << norm_ids[i] << " " << self_coop_levels[i] << " ";
    for (size_t j = 0; j < p_fix.Cols(); j++) {
      out << p_fix(i,j) << " ";
    }
    out << std::endl;
  }
}

struct Input {
  EvolPrivRepGame::Parameters evoparam;
  std::vector<std::pair<double, double>> benefit_sigma_in_vec;

  Input() : evoparam(), benefit_sigma_in_vec() {};

  NLOHMANN_DEFINE_TYPE_INTRUSIVE_WITH_DEFAULT(Input, evoparam, benefit_sigma_in_vec);
};


int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);

  int my_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  using namespace nlohmann;
  enum NormSet {
    LeadingEight,
    SecondOrder,
    ThirdOrder,
    DualSecondOrder
  };

  json j = json::object();
  bool debug_mode = false;
  NormSet norm_set = LeadingEight;
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
    else if (std::string(argv[i]) == "-d") {
      debug_mode = true;
    }
    else if (std::string(argv[i]) == "--norm-set") {
      std::string set_name = argv[++i];
      if (set_name == "leading_eight") {
        norm_set = LeadingEight;
      }
      else if (set_name == "second_order") {
        norm_set = SecondOrder;
      }
      else if (set_name == "third_order") {
        norm_set = ThirdOrder;
      }
      else if (set_name == "dual_second_order") {
        norm_set = DualSecondOrder;
      }
      else {
        std::cerr << "unknown norm set: " << set_name << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }
    else {
      std::cerr << "unknown option: " << argv[i] << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  Input input = j.get<Input>();

  if (my_rank == 0) {
    std::cerr << "params: " << nlohmann::json(input) << std::endl;
    std::cerr << "norm_set: " << norm_set << std::endl;
  }

  if (debug_mode) {
    // debug mode
    if (my_rank == 0) {
      std::cerr << "debug mode" << std::endl;
    }
    MPI_Finalize();
    return 0;
  }

  // measure elapsed time
  auto start = std::chrono::system_clock::now();

  std::vector<Norm> norms;
  std::string norm_set_name;
  if (norm_set == LeadingEight) {
    norms = {Norm::L1(), Norm::L2(), Norm::L3(), Norm::L4(), Norm::L5(), Norm::L6(), Norm::L7(), Norm::L8(), Norm::AllC(), Norm::AllD()};
    norm_set_name = "leading_eight";
  }
  else if (norm_set == SecondOrder) {
    norms = Norm::Deterministic2ndOrderWithoutR2Norms();
    norm_set_name = "second_order";
  }
  else if (norm_set == ThirdOrder) {
    norms = Norm::Deterministic3rdOrderWithoutR2Norms();
    norm_set_name = "third_order";
  }
  else if (norm_set == DualSecondOrder) {
    norms = Norm::Deterministic2ndOrderWithR2Norms();
    norm_set_name = "dual_second_order";
  }

  std::vector<double> self_coop_levels;
  p_fix_vec_t p_fix_vec;
  CalculateFixationProbs(input.evoparam, input.benefit_sigma_in_vec, norms, self_coop_levels, p_fix_vec);

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
      auto [benefit, sigma_in] = input.benefit_sigma_in_vec.at(n);
      WriteInMsgpack(filename, input.evoparam, benefit, sigma_in, norm_set_name, norm_ids, self_coop_levels, p_fix_vec[n]);
      // write txt file as well
      std::ofstream fout("fixation_probs_" + std::to_string(n) + ".dat");
      PrintFixationProbsInText(fout, norm_set_name, norm_ids, self_coop_levels, p_fix_vec[n]);
    }
  }

  MPI_Finalize();
  return 0;
}

