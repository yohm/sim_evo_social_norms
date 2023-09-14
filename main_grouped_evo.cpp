#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <chrono>
#include <random>
#include <nlohmann/json.hpp>
#include <icecream.hpp>
#include "Vector2d.hpp"
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

class GroupedEvoGame {
public:
  class Parameters {
    public:
    Parameters() : M(100), T_init(1e6), T_measure(1e7), seed(123456789ull), benefit(5.0), sigma_out(5.0), mut_r(0.1) {}
    size_t M;   // number of groups
    size_t T_init;
    size_t T_measure;
    uint64_t seed;
    double benefit;
    double sigma_out;
    double mut_r;   // relative mutation rate

    NLOHMANN_DEFINE_TYPE_INTRUSIVE_WITH_DEFAULT(GroupedEvoGame::Parameters, M, T_init, T_measure, seed, benefit, sigma_out, mut_r)
  };

  explicit GroupedEvoGame(const Parameters& _prm, const std::vector<Norm>& _norms) :
    prm(_prm), norms(_norms)
  {
    if (_norms.empty()) {
      throw std::runtime_error("norms is empty");
    }
    rng.seed(prm.seed);
    std::uniform_int_distribution<size_t> uni_int(0, norms.size()-1);
    species.reserve(prm.M);
    for (size_t i = 0; i < prm.M; i++) {
      species.emplace_back( uni_int(rng) );
    }
    histogram.assign(norms.size(), 0);
  }
  Parameters prm;
  std::vector<size_t> species;   //  species[i] : species index at group i
  std::vector<Norm> norms;
  std::mt19937_64 rng;
  Vector2d<double> fixation_prob_cache;
  std::vector<double> self_coop_level_cache;
  std::vector<size_t> histogram;  // histogram of species during evolution

  void SetFixationProbsCache(const Vector2d<double>& fixation_probs) {
    fixation_prob_cache = fixation_probs;
  }
  void SetSelfCoopLevelCache(const std::vector<double>& self_c_probs) {
    self_coop_level_cache = self_c_probs;
  }
  void Update() {
    // one Monte Carlo sweep
    for (int t = 0; t < prm.M; t++) {
      std::uniform_int_distribution<size_t> d0(0, prm.M-1);
      size_t res_index = d0(rng);
      UpdateGroup(res_index);
    }
  }
private:
  void UpdateGroup(size_t res_index) {
    // focal species : species[g]
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    size_t resident = species[res_index];

    if (uni(rng) < prm.mut_r) {
      // mutation
      std::uniform_int_distribution<size_t> d1(1, norms.size()-1);
      size_t mutant = (resident + d1(rng)) % norms.size();
      double f = IntraGroupFixationProb(mutant, resident);
      if (uni(rng) < f) {
        species[res_index] = mutant;
      }
    }
    else {
      std::uniform_int_distribution<size_t> d1(1, prm.M-1);
      size_t mig_index = static_cast<size_t>(res_index + d1(rng)) % prm.M;
      size_t immigrant = species[mig_index];
      double p = InterGroupImitationProb(immigrant, resident);
      double f = IntraGroupFixationProb(immigrant, resident);
      double prob = p*f;
      if (uni(rng) < prob) {
        species[res_index] = immigrant;
      }
    }
  }
  double IntraGroupFixationProb(size_t mutant, size_t resident) const {
    return fixation_prob_cache(resident, mutant);
  }
  double InterGroupImitationProb(size_t immigrant, size_t resident) const {
    double benefit = prm.benefit;
    double pi_resident  = (benefit - 1.0) * self_coop_level_cache[resident];
    double pi_immigrant = (benefit - 1.0) * self_coop_level_cache[immigrant];
    // f_{A\to B} = { 1 + \exp[ \sigma_out (s_A - s_B) ] }^{-1}
    return 1.0 / (1.0 + std::exp(prm.sigma_out * (pi_resident - pi_immigrant) ));
  }
public:
  size_t NumGroupsOf(const Norm& n) const {
    size_t count = 0;
    for (size_t s : species) {
      if (norms[s] == n) {
        count++;
      }
    }
    return count;
  }
  void UpdateHistogram() {
    for (unsigned long s : species) {
      histogram[s] += 1;
    }
  }
  void ResetHistogram() {
    std::fill(histogram.begin(), histogram.end(), 0);
  }
  [[nodiscard]] std::vector<double> NormalizedHistogram() const {
    std::vector<double> normed(histogram.size(), 0.0);
    double sum = 0.0;
    for (unsigned long h : histogram) {
      sum += static_cast<double>(h);
    }
    for (size_t i = 0; i < histogram.size(); i++) {
      normed[i] = static_cast<double>(histogram[i]) / sum;
    }
    return normed;
  }
  [[nodiscard]] double AverageCoopLevelOverHistogram() const {
    std::vector<double> n_histo = NormalizedHistogram();
    double sum = 0.0;
    for (size_t i = 0; i < n_histo.size(); i++) {
      sum += n_histo[i] * self_coop_level_cache[i];
    }
    return sum;
  }
  [[nodiscard]] double CurrentCooperationLevel() const {
    double sum = 0.0;
    for (const auto& s: species) {
      sum += self_coop_level_cache[s];
    }
    return sum / static_cast<double>(species.size());
  }
};

void PrintProgress(double progress) {
  int barWidth = 70;
  std::cerr << "[";
  int pos = static_cast<int>(barWidth * progress);
  for (int i = 0; i < barWidth; ++i) {
    if (i < pos) std::cerr << "=";
    else if (i == pos) std::cerr << ">";
    else std::cerr << " ";
  }
  std::cerr << "] " << int(progress * 100.0) << " %\r";
  std::cerr.flush();
  if (progress >= 1.0) {
    std::cerr << std::endl;
    return;
  }
}

int main(int argc, char* argv[]) {
  // run evolutionary simulation in group-structured population

  using json = nlohmann::json;

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
      args.emplace_back(argv[i]);  // for backward compatibility (do nothing
    }
  }
  GroupedEvoGame::Parameters params = j.get<GroupedEvoGame::Parameters>();
  std::cerr << "parameters:" << std::endl;
  std::cerr << nlohmann::json(params) << std::endl;

  if (args.size() != 1) {
    std::cerr << "Error: invalid number of arguments" << std::endl;
    std::cerr << "Usage: " << argv[0] << " <input_msgpack_file> -j [parameter.json]" << std::endl;
    return 1;
  }


  // load fixation probabilities from the input file
  json j_in = LoadMsgpackFile(argv[1]);
  std::cerr << "json loaded" << std::endl;

  std::vector<Norm> norms;
  for (auto& norm_id_j : j_in["norm_ids"]) {
    int id = norm_id_j.get<int>();
    Norm n = Norm::ConstructFromID(id);
    norms.emplace_back(n);
  }
  std::cerr << "  # of norms: " << norms.size() << std::endl;

  const size_t N_NORMS = norms.size();
  Vector2d<double> p_fix(N_NORMS, N_NORMS, 0.0);
  std::vector<double> self_coop_levels(N_NORMS, 0.0);
  for (size_t i = 0; i < p_fix.size(); ++i) {
    p_fix._data[i] = j_in["p_fix"][i].get<double>();
  }
  for (size_t i = 0; i < self_coop_levels.size(); ++i) {
    self_coop_levels[i] = j_in["self_coop_levels"][i].get<double>();
  }
  double benefit = j_in["benefit"].get<double>();
  // check consistency of benefit
  if ( std::abs(benefit - params.benefit) > 1.0e-8 ) {
    std::cerr << "Error: benefit in the input file is inconsistent with the parameter file" << std::endl;
    std::cerr << j_in["benefit"] << std::endl;
    return 1;
  }
  // IC(p_fix._data, self_coop_levels);

  GroupedEvoGame evo(params, norms);
  evo.SetFixationProbsCache(p_fix);
  evo.SetSelfCoopLevelCache(self_coop_levels);

  std::ofstream tout("timeseries.dat");
  std::vector<Norm> norms_to_measure = {Norm::L1(), Norm::ConstructFromID(765130), Norm::L3(), Norm::ConstructFromID(769227)};
  auto plot_time_series = [&tout,&norms_to_measure](size_t t, const GroupedEvoGame& evo) {
    tout << t << ' ' << evo.CurrentCooperationLevel();
    for (const auto& n: norms_to_measure) {
      tout << ' ' << static_cast<double>(evo.NumGroupsOf(n)) / static_cast<double>(evo.prm.M);
    }
    tout << std::endl;
  };
  size_t print_interval = (params.T_init + params.T_measure) / 1000;
  size_t progress_interval = (params.T_init + params.T_measure) / 100;


  for (size_t t = 0; t < params.T_init; t++) {
    if (t % progress_interval == 0) {
      PrintProgress(static_cast<double>(t) / static_cast<double>(params.T_init+params.T_measure));
    }
    evo.Update();
    if (t % print_interval == 0) {
      plot_time_series(t, evo);
    }
  }
  for (size_t t = params.T_init; t < params.T_measure + params.T_init; t++) {
    if (t % progress_interval == 0) {
      PrintProgress(static_cast<double>(t) / static_cast<double>(params.T_init+params.T_measure));
    }
    evo.Update();
    evo.UpdateHistogram();
    if (t % print_interval == 0) {
      plot_time_series(t, evo);
    }
  }
  PrintProgress(1.0);
  tout.close();

  auto n_histo = evo.NormalizedHistogram();
  {
    double c_level = evo.AverageCoopLevelOverHistogram();
    std::cerr << "overall_c_prob: " << c_level << std::endl;
    std::ofstream jout("_output.json");
    nlohmann::json j;
    j["cooperation_level"] = evo.AverageCoopLevelOverHistogram();
    for (const Norm& n: norms_to_measure) {
      // get an index of n in norms
      auto it = std::find(norms.begin(), norms.end(), n);
      if (it != norms.end()) {
        size_t idx = std::distance(norms.begin(), it);
        j["norms"][std::to_string(n.ID())] = n_histo[idx];
      }
      else {
        std::cerr << "Error: norm " << n.ID() << " not found" << std::endl;
      }
    }

    jout << j;
    jout.close();
  }

  std::vector<Norm> l8 = {Norm::L1(), Norm::L2(), Norm::L3(), Norm::L4(), Norm::L5(), Norm::L6(), Norm::L7(), Norm::L8()};
  std::vector<Norm> l8v;
  for (Norm l : l8) {
    constexpr Reputation B = Reputation::B;
    l.P.SetCProb(B, B, 1.0 - l.P.CProb(B, B));
    l8v.push_back(l);
  }

  auto l8_variant_name = [&l8v](const Norm& n) -> std::string {
    // find n from l8v, return its index
    auto it = std::find(l8v.begin(), l8v.end(), n);
    if (it != l8v.end()) {
      size_t idx = std::distance(l8v.begin(), it) + 1;
      return "L" + std::to_string(idx) + "v";
    }
    return "";
  };

  // sort histo in descending order with its index
  std::vector<std::pair<double, size_t>> histo_sorted;
  for (size_t i = 0; i < n_histo.size(); ++i) {
    histo_sorted.emplace_back(n_histo[i], i);
  }
  std::sort(histo_sorted.begin(), histo_sorted.end(), std::greater<>());
  std::ofstream histo_out("histo_norms.dat");
  for (size_t i = 0; i < histo_sorted.size(); ++i) {
    size_t idx = histo_sorted[i].second;
    int nid = norms[idx].ID();
    std::string type = "other";
    if (norms[idx].P.ID() == 0) { type = "AllD"; }
    else if (norms[idx].P.ID() == 15) { type = "AllC"; }
    else if (!norms[idx].GetName().empty()) { type = norms[idx].GetName(); }
    else if (!l8_variant_name(norms[idx]).empty()) {
      type = l8_variant_name(norms[idx]);
    }
    histo_out << idx << ' ' << nid << ' ' << histo_sorted[i].first << ' ' << self_coop_levels[idx] << ' ' << type << std::endl;
    if (histo_sorted[i].first < 0.01 && i > 8) {
      break;
    }
  }
  histo_out.close();

  return 0;
}
