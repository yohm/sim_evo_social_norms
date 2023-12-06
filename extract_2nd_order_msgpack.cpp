#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <chrono>
#include <random>
#include <nlohmann/json.hpp>
#include "Norm.hpp"


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
  // run evolutionary simulation in homogeneous population

  using json = nlohmann::json;

  std::string msgpack_path, out_path;

  if (argc == 3) {
    msgpack_path = argv[1];
    out_path = argv[2];
  }
  else {
    std::cerr << "Usage: " << argv[0] << " <input_msgpack_file> <out_msgpack_path>" << std::endl;
    return 1;
  }

  // load fixation probabilities from the input file
  json j_in = LoadMsgpackFile(msgpack_path);
  std::cerr << "json loaded" << std::endl;

  json j_out = j_in;

  std::vector<size_t> indices;
  std::vector<int> new_norm_ids;

  for (size_t i = 0; i < j_in["norm_ids"].size(); i++) {
    int id = j_in["norm_ids"][i].get<int>();
    Norm n = Norm::ConstructFromID(id);
    if (n.IsSecondOrder()) {
      indices.emplace_back(i);
      new_norm_ids.emplace_back(id);
    }
  }
  std::cerr << "  # of norms: " << indices.size() << std::endl;

  const size_t N_NORMS = j_in["norm_ids"].size();

  std::vector<double> p_fix_flattened;
  std::vector<double> self_coop_levels;

  for (size_t i = 0; i < N_NORMS; ++i) {
    // if i is included in indices
    if (std::find(indices.begin(), indices.end(), i) != indices.end()) {
      for (size_t j = 0; j < N_NORMS; ++j) {
        // if j is included in indices
        if (std::find(indices.begin(), indices.end(), j) != indices.end()) {
          p_fix_flattened.push_back( j_in["p_fix"][i * N_NORMS + j].get<double>() );
        }
      }
      self_coop_levels.push_back( j_in["self_coop_levels"][i].get<double>() );
    }
  }

  j_out["norm_set"] = "second_order";
  j_out["norm_ids"] = new_norm_ids;
  j_out["p_fix"] = p_fix_flattened;
  j_out["self_coop_levels"] = self_coop_levels;

  // write file
  std::ofstream ofs(out_path, std::ios::binary);
  std::vector<std::uint8_t> v = nlohmann::json::to_msgpack(j_out);
  ofs.write(reinterpret_cast<char*>(v.data()), v.size());

  return 0;
}
