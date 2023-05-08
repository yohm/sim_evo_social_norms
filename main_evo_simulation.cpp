#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <chrono>
#include <random>
#include <nlohmann/json.hpp>
#include "Norm.hpp"
#include "PrivRepGame.hpp"


constexpr Reputation B = Reputation::B, G = Reputation::G;
constexpr Action C = Action::C, D = Action::D;

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

void SimulateWellMixedPopulation(const std::vector<Norm>& norms, const vector2d<double>& p_fix,
                                 const std::vector<double>& self_coop_levels, uint64_t seed,
                                 size_t T_init, size_t T_measure) {

  // initialize random number generator
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> uni(0.0, 1.0);
  auto r01 = [&uni, &rng] { return uni(rng); };

  // initialize population
  size_t resident = r01() * norms.size();

  for (size_t t = 0; t < T_init; t++) {
    size_t mut = (resident + 1 + static_cast<size_t>(r01() * (norms.size()-1))) % norms.size();
    assert(resident != mut);
    if (r01() < p_fix(resident, mut)) {
      resident = mut;
    }

    constexpr size_t interval = 100ul;
    if (t % interval == 0) {
      std::cerr << t << ' ' << resident << ' ' << self_coop_levels[resident] << std::endl;
    }
  }
}

int main(int argc, char* argv[]) {
  // run evolutionary simulation in group-structured population

  using namespace nlohmann;

  size_t T_init = 1e4;
  size_t T_measure = 1e4;
  uint64_t seed = 123456789ull;

  // load fixation probabilities from the input file
  json j_in = LoadMsgpackFile("fixation_probs_3species.msgpack");
  std::cerr << j_in << std::endl;

  std::vector<Norm> norms;
  for (auto& norm_id_j : j_in["norm_ids"]) {
    int id = norm_id_j.get<int>();
    Norm n = Norm::ConstructFromID(id);
    norms.emplace_back(n);
  }

  const size_t N_NORMS = norms.size();
  vector2d<double> p_fix(N_NORMS, N_NORMS, 0.0);
  std::vector<double> self_coop_levels(N_NORMS, 0.0);
  std::cerr << j_in["p_fix"][0].get<double>() << std::endl;
  for (size_t i = 0; i < p_fix.size(); ++i) {
    p_fix._data[i] = j_in["p_fix"][i].get<double>();
  }
  for (size_t i = 0; i < self_coop_levels.size(); ++i) {
    self_coop_levels[i] = j_in["self_coop_levels"][i].get<double>();
  }
  IC(p_fix._data, self_coop_levels);

  SimulateWellMixedPopulation(norms, p_fix, self_coop_levels, seed, T_init, T_measure);

  return 0;
}

