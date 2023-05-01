#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <chrono>
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
  // bytes.reserve(length);
  bytes.resize(length);

  // read bytes into vector
  char* buffer = new char[length];
  fin.read(bytes.data(), length);

    // unpack msgpack
  nlohmann::json j = nlohmann::json::from_msgpack(bytes);

  std::cerr << j;
  return j;
}


int main(int argc, char* argv[]) {
  // run evolutionary simulation in group-structured population

  using namespace nlohmann;

  size_t M = 1; // number of groups
  size_t t_init = 1e4;
  size_t t_measure = 1e4;
  uint64_t seed = 123456789ull;

  // load fixation probabilities from the input file
  json j_in = LoadMsgpackFile("fixation_probs.msgpack");

  // calculate intra-group fixation probabilities
  // constexpr size_t N_NORMS = 4096;
  // vector2d<double> p_fix(N_NORMS, N_NORMS, 0.0);
  // std::vector<double> self_coop_levels(N_NORMS, 0.0);

  return 0;
}

