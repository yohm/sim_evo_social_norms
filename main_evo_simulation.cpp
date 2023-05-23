#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <chrono>
#include <random>
#include <nlohmann/json.hpp>
#include "Vector2d.hpp"
#include "Norm.hpp"
#include "PrivRepGame.hpp"


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
    Parameters() : M(100), T_init(1e6), T_measure(1e7), seed(123456789ull), benefit(5.0), sigma_out(1.0), mut_r(0.1) {}
    size_t M;   // number of groups
    size_t T_init;
    size_t T_measure;
    uint64_t seed;
    double benefit;
    double sigma_out;
    double mut_r;   // relative mutation rate
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
    histogram.resize(norms.size(), 0);
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
  double IntraGroupFixationProb(size_t mutant, size_t resident) {
    return fixation_prob_cache(resident, mutant);
  }
  double InterGroupImitationProb(size_t immigrant, size_t resident) {
    double pi_resident  = (prm.benefit - 1.0) * self_coop_level_cache[resident];
    double pi_immigrant = (prm.benefit - 1.0) * self_coop_level_cache[immigrant];
    // f_{A\to B} = { 1 + \exp[ \sigma_out (s_A - s_B) ] }^{-1}
    return 1.0 / (1.0 + std::exp(prm.sigma_out * (pi_resident - pi_immigrant) ));
  }
  void TakeHistogram() {
    for (size_t i = 0; i < species.size(); i++) {
      histogram[species[i]] += 1;
    }
  }
  void ResetHistogram() {
    std::fill(histogram.begin(), histogram.end(), 0);
  }
  std::vector<double> NormalizedHistogram() const {
    std::vector<double> normed(histogram.size(), 0.0);
    double sum = 0.0;
    for (size_t i = 0; i < histogram.size(); i++) {
      sum += histogram[i];
    }
    for (size_t i = 0; i < histogram.size(); i++) {
      normed[i] = histogram[i] / sum;
    }
    return normed;
  }
  double AverageCoopLevel() const {
    std::vector<double> n_histo = NormalizedHistogram();
    double sum = 0.0;
    for (size_t i = 0; i < n_histo.size(); i++) {
      sum += n_histo[i] * self_coop_level_cache[i];
    }
    return sum;
  }
};

int main(int argc, char* argv[]) {
  // run evolutionary simulation in group-structured population

  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <input_msgpack_file>" << std::endl;
    exit(1);
  }

  using namespace nlohmann;

  // load fixation probabilities from the input file
  json j_in = LoadMsgpackFile(argv[1]);
  std::cerr << "json loaded" << std::endl;

  std::vector<Norm> norms;
  for (auto& norm_id_j : j_in["norm_ids"]) {
    int id = norm_id_j.get<int>();
    Norm n = Norm::ConstructFromID(id);
    norms.emplace_back(n);
  }

  const size_t N_NORMS = norms.size();
  Vector2d<double> p_fix(N_NORMS, N_NORMS, 0.0);
  std::vector<double> self_coop_levels(N_NORMS, 0.0);
  for (size_t i = 0; i < p_fix.size(); ++i) {
    p_fix._data[i] = j_in["p_fix"][i].get<double>();
  }
  for (size_t i = 0; i < self_coop_levels.size(); ++i) {
    self_coop_levels[i] = j_in["self_coop_levels"][i].get<double>();
  }
  // IC(p_fix._data, self_coop_levels);

  // SimulateWellMixedPopulation(norms, p_fix, self_coop_levels, seed, T_init, T_measure);

  GroupedEvoGame::Parameters params;
  GroupedEvoGame evo(params, norms);
  evo.SetFixationProbsCache(p_fix);
  evo.SetSelfCoopLevelCache(self_coop_levels);
  for (size_t t = 0; t < params.T_init; t++) {
    evo.Update();
  }
  for (size_t t = 0; t < params.T_measure; t++) {
    size_t interval = params.T_measure / 100;
    if (t % interval == 0) {
      std::cerr << "t = " << t << std::endl;
    }
    evo.Update();
    evo.TakeHistogram();
  }
  auto n_histo = evo.NormalizedHistogram();
  // Print n_histo
  for (size_t i = 0; i < n_histo.size(); ++i) {
    std::cout << i << ' ' << norms[i].ID() << ' ' << n_histo[i] << ' ' << self_coop_levels[i] << std::endl;
  }
  // print average cooperation level
  std::cout << "overall_c_prob: " << evo.AverageCoopLevel() << std::endl;
  std::cout << "histo: " << std::endl;

  // sort histo in descending order with its index
  std::vector<std::pair<double, size_t>> histo_sorted;
  for (size_t i = 0; i < n_histo.size(); ++i) {
    histo_sorted.emplace_back(std::make_pair(n_histo[i], i));
  }
  std::sort(histo_sorted.begin(), histo_sorted.end(), std::greater<std::pair<double, size_t>>());
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

  return 0;
}

