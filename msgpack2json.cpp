#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>


int main(int argc, char* argv[]) {

  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <input_msgpack_file>" << std::endl;
    return 1;
  }

  std::ifstream fin(argv[1], std::ios::binary);
  std::vector<char> bytes;
  // read the fin into the vector
  if (!fin) {
    std::cerr << "Error opening file " << argv[1] << std::endl;
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

  // print json with indentation
  std::cout << j.dump(2) << std::endl;

  return 0;
}