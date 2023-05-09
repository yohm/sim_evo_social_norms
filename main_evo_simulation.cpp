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
  std::uniform_int_distribution<size_t> uni_int(0, norms.size()-1);
  auto r01 = [&uni, &rng] { return uni(rng); };

  // initialize population
  size_t resident = r01() * norms.size();

  size_t count = 0;
  double overall_c_prob = 0.0;
  std::vector<double> histo(norms.size(), 0.0);

  for (size_t t = 0; t < T_init+T_measure; t++) {
    size_t mut = (resident + 1 + uni_int(rng)) % norms.size();
    assert(resident != mut);
    if (r01() < p_fix(resident, mut)) {
      resident = mut;
    }

    const size_t interval = T_measure / 100;
    if (t % interval == 0) {
      std::cerr << t << ' ' << resident << ' ' << self_coop_levels[resident] << std::endl;
    }

    if (t > T_init) {
      overall_c_prob += self_coop_levels[resident];
      histo[resident] += 1.0;
      count++;
    }
  }

  // normalize cooperation_level & histogram
  overall_c_prob /= count;
  for (size_t i = 0; i < histo.size(); ++i) {
      histo[i] /= count;
  }
  // sort histo in descending order with its index
  std::vector<std::pair<double, size_t>> histo_sorted;
  for (size_t i = 0; i < histo.size(); ++i) {
      histo_sorted.emplace_back(std::make_pair(histo[i], i));
  }
  std::sort(histo_sorted.begin(), histo_sorted.end(), std::greater<std::pair<double, size_t>>());

  std::cout << "overall_c_prob: " << overall_c_prob << std::endl;
  std::cout << "histo: " << std::endl;
  for (size_t i = 0; i < histo_sorted.size(); ++i) {
    size_t idx = histo_sorted[i].second;
    int nid = norms[idx].ID();
    std::string type = "other";
    if (norms[idx].P.ID() == 0) { type = "AllD"; }
    else if (norms[idx].P.ID() == 15) { type = "AllC"; }
    else if (!norms[idx].GetName().empty()) { type = norms[idx].GetName(); }
    std::cout << idx << ' ' << nid << ' ' << histo_sorted[i].first << ' ' << self_coop_levels[idx] << ' ' << type << std::endl;
    if (histo_sorted[i].first < 0.01 && i > 20) {
      break;
    }
  }
}

int main(int argc, char* argv[]) {
  // run evolutionary simulation in group-structured population

  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <input_msgpack_file>" << std::endl;
    exit(1);
  }

  using namespace nlohmann;

  size_t T_init = 1e6;
  size_t T_measure = 1e8;
  uint64_t seed = 123456789ull;

  // load fixation probabilities from the input file
  json j_in = LoadMsgpackFile(argv[1]);
  // std::cerr << j_in << std::endl;

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
  // IC(p_fix._data, self_coop_levels);

  SimulateWellMixedPopulation(norms, p_fix, self_coop_levels, seed, T_init, T_measure);

  return 0;
}

