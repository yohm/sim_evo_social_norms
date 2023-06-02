#include <iostream>
#include <fstream>
#include <chrono>
#include <regex>
#include <icecream.hpp>
#include <nlohmann/json.hpp>
#include "PrivRepGame.hpp"
#include "Parameters.hpp"


using timeseries_t = std::vector<std::pair<int,double>>;

timeseries_t TimeSeries(PrivateRepGame& prg, const Parameters& params) {
  size_t total_steps = params.n_init + params.n_steps;
  size_t num_prints = 100;
  size_t interval = total_steps / num_prints;

  std::vector< std::pair<int,double>> timeseries;
  for (size_t i = 0; i < num_prints; i++) {
    prg.Update(interval, params.q, params.mu_percept, false);
    double pc = prg.NormCooperationLevels()[0][0];
    timeseries.emplace_back((i+1)*interval, pc);
    prg.ResetCounts();
  }

  return timeseries;
}

std::vector<timeseries_t> PolymorphicTimeSeries(const Norm& n1, const Norm& n2, const Parameters& params) {
  PrivateRepGame prg1({{n1, params.N}}, params.seed);
  auto t1 = TimeSeries(prg1, params);

  PrivateRepGame prg2({{n1, params.N/2}, {n2, params.N/2}}, params.seed);
  auto t2 = TimeSeries(prg2, params);

  PrivateRepGame prg3({{n2, params.N}}, params.seed);
  auto t3 = TimeSeries(prg3, params);

  return {t1, t2, t3};
}

int main(int argc, char *argv[]) {

  std::vector<std::string> args;
  bool swap_gb = false;
  nlohmann::json j = nlohmann::json::object();
  // -j param.json : set parameters used for evolutionary simulation by json file
  // -s : swap good and bad
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
    else if (std::string(argv[i]) == "-s") {
      swap_gb = true;
    }
    else {
      args.emplace_back(argv[i]);
    }
  }

  auto start = std::chrono::high_resolution_clock::now();

  if (args.size() == 1) {
    Norm n = Norm::ParseNormString(args.at(0), swap_gb);
    std::cerr << n.Inspect();
    Parameters params = j.get<Parameters>();
    std::cerr << nlohmann::json(params).dump(2) << std::endl;
    PrivateRepGame prg({{n, params.N}}, params.seed);
    auto t_series = TimeSeries(prg, params);
    for (auto& t : t_series) {
      std::cout << t.first << " " << t.second << std::endl;
    }
  }
  else if (args.size() == 2) {  // if two arguments are given, direct competition between two norms are shown
    Norm n1 = Norm::ParseNormString(args.at(0), swap_gb);
    Norm n2 = Norm::ParseNormString(args.at(1), swap_gb);
    Parameters params = j.get<Parameters>();
    std::cerr << nlohmann::json(params).dump(2) << std::endl;
    auto vt = PolymorphicTimeSeries(n1, n2, params);
    for (size_t t = 0; t < vt[0].size(); t++) {
      std::cout << vt[0][t].first << " " << vt[0][t].second << " " << vt[1][t].second << " " << vt[2][t].second << std::endl;
    }
  }
  else {
    std::cerr << "Usage: " << argv[0] << " [-j param.json] [-s] norm_string" << std::endl;
    std::cerr << "       " << argv[0] << " [-j param.json] [-s] norm_string1 norm_string2" << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << "  -j param.json : set parameters used for evolutionary simulation by json file" << std::endl;
    std::cerr << "  -s : swap good and bad" << std::endl;
    std::cerr << "  norm_string : string representation of a norm" << std::endl;
    std::cerr << "Default parameters:" << std::endl;
    std::cerr << "  " << nlohmann::json(Parameters()).dump(2) << std::endl;
    return 1;
  }

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cerr << "Elapsed time: " << elapsed.count() << " s\n";

  return 0;
}