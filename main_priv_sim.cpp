#include <iostream>
#include <cassert>
#include <icecream.hpp>
#include <vector>
#include <random>
#include <mpi.h>
#include <omp.h>
#include <nlohmann/json.hpp>
#include <caravan.hpp>
#include "Norm.hpp"
#include "PrivRepGame.hpp"
#include "Parameters.hpp"


constexpr Reputation B = Reputation::B, G = Reputation::G;
constexpr Action C = Action::C, D = Action::D;

double SelfCoopLevel(const Norm& norm, const Parameters& params) {
  PrivateRepGame prg({{norm, params.N}}, params.seed);
  prg.Update(params.n_init, params.q, params.mu_percept, false);
  prg.ResetCounts();
  prg.Update(params.n_steps, params.q, params.mu_percept, true);
  return prg.NormCooperationLevels()[0][0];
}

const EvolPrivRepGameAllCAllD* p_evol = nullptr;
void InitializeEvol(const Parameters& params) {
  EvolPrivRepGame::SimulationParameters evoparams = params.ToEvolParams();
  p_evol = new EvolPrivRepGameAllCAllD(evoparams, params.benefit, params.beta);
}

const EvolPrivRepGameAllCAllD GetEvol() {
  if (!p_evol) {
    throw std::runtime_error("EvolPrivRepGameAllCAllD is not initialized!");
  }
  return *p_evol;
}

// return <ID, cooperation level> for a given norm
double EqCooperationLevel(const Norm& norm) {
  const auto evol = GetEvol();
  auto pc_rho_eq = evol.EquilibriumCoopLevelAllCAllD(norm);
  double self_cooperation_level = std::get<0>(pc_rho_eq);
  auto eq = std::get<2>(pc_rho_eq);
  return self_cooperation_level * eq[0] + 1.0 * eq[1];
}

std::vector<Norm> LocalMutants(const Norm& norm) {
  std::vector<Norm> local_mutants;

  for (int i = 0; i < 20; i++) {
    auto serialized = norm.Serialize();
    const double delta = 1.0;
    if (serialized[i] <= 0.5) {
      serialized[i] += delta;
      if (serialized[i] > 1.0) {
        serialized[i] = 1.0;
      }
    }
    else {
      serialized[i] -= delta;
      if (serialized[i] < 0.0) {
        serialized[i] = 0.0;
      }
    }
    Norm mutant = Norm::FromSerialized(serialized);
    local_mutants.push_back(mutant);
  }

  return local_mutants;
}

std::pair<double,double> EqCooperationLevelWithLocalMutants(const Norm& norm, const Parameters& params) {
  EvolPrivRepGame::SimulationParameters evoparams = params.ToEvolParams();
  std::vector<Norm> norms = {norm, Norm::AllC(), Norm::AllD()};
  auto local_mutants = LocalMutants(norm);
  norms.insert(norms.end(), local_mutants.begin(), local_mutants.end());
  EvolPrivRepGame evol(evoparams);
  auto fixation_probs = evol.FixationProbabilities(norms, params.benefit, params.beta);
  std::vector<double> eq = EvolPrivRepGame::EquilibriumPopulationLowMut(fixation_probs);
  std::vector<double> self_coop_levels;
  for (const auto& norm : norms) {
    self_coop_levels.push_back(evol.SelfCooperationLevel(norm));
  }

  double eq_coop_level = 0.0;
  for (size_t i = 0; i < norms.size(); i++) {
    eq_coop_level += self_coop_levels[i] * eq[i];
  }
  return std::make_pair(eq_coop_level, eq[0]);
}

std::pair<double,Norm> MostRiskyMutant(const Norm& norm, const Parameters& params) {
  EvolPrivRepGame::SimulationParameters evoparams = params.ToEvolParams();

  double min_eq_population = 1.0;
  Norm most_risky_mutant = Norm::AllC();
  std::vector<Norm> mutants = {Norm::AllC(), Norm::AllD()};
  auto local_mutants = LocalMutants(norm);
  mutants.insert(mutants.end(), local_mutants.begin(), local_mutants.end());
  for (const auto& mutant : mutants) {
    EvolPrivRepGame evol(evoparams);
    auto fixation_probs = evol.FixationProbabilities({norm, mutant}, params.benefit, params.beta);
    std::vector<double> eq = EvolPrivRepGame::EquilibriumPopulationLowMut(fixation_probs);
    if (eq[0] < min_eq_population) {
      min_eq_population = eq[0];
      most_risky_mutant = mutant;
    }
    if (min_eq_population < 1.0e-2) {
      break;
    }
  }

  return std::make_pair(min_eq_population, most_risky_mutant);
}


int main(int argc, char** argv) {
  using namespace nlohmann;

  MPI_Init(&argc, &argv);
  int my_rank = 0, num_procs = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  std::vector<std::string> args;
  json j = json::object();
  bool without_R2 = false;
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
    else if (std::string(argv[i]) == "--without-R2") {
      without_R2 = true;
    }
    else {
      std::cerr << "unknown option: " << argv[i] << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }
  Parameters params = j.get<Parameters>();
  InitializeEvol(params);

  if (my_rank == 0) {
    std::cerr << "params: " << json(params).dump(2) << "\n";
    std::cerr << "without_R2: " << without_R2 << "\n";
    std::cerr << "num_procs: " << num_procs << "\n";
    std::cerr << "num_threads: " << omp_get_max_threads() << "\n";
  }

  std::function<void(caravan::Queue&)> on_init = [without_R2](caravan::Queue& q) {
    int i_max = 256, j_min = 0, j_max = 256;
    if (without_R2) {
      j_min = AssessmentRule::KeepRecipient().ID();
      j_max = j_min + 1;
    }
    for (int i = 0; i < i_max; i++) {
      for (int j = j_min; j < j_max; j++) {
        json input = {i, j};
        q.Push(input);
      }
    }
  };

  std::vector< std::tuple<int,double,double> > results;
  std::function<void(int64_t, const json&, const json&, caravan::Queue&)> on_result_receive = [&results](int64_t task_id, const json& input, const json& output, caravan::Queue& q) {
    if (task_id % 1000 == 0) {
      std::cerr << "task: " << task_id << " has finished: input: " << input << ", output: " << output << "\n";
    }
    for (auto result: output) {
      results.emplace_back(result.at(0).get<int>(), result.at(1).get<double>(), result.at(2).get<double>());
    }
  };

  std::function<json(const json& input)> do_task = [&params](const json& input) {
    json output = json::array();
    int R1_id = input.at(0).get<int>();
    int R2_id = input.at(1).get<int>();
    AssessmentRule R1 = AssessmentRule::MakeDeterministicRule(R1_id);
    AssessmentRule R2 = AssessmentRule::MakeDeterministicRule(R2_id);
    for (int i = 0; i < 16; i++) {
      ActionRule P = ActionRule::MakeDeterministicRule(i);
      Norm norm(R1, R2, P);
      if (norm.ID() < norm.SwapGB().ID()) {
        continue;
      }

      double self_coop_level = SelfCoopLevel(norm, params);
      auto p = MostRiskyMutant(norm, params);
      output.emplace_back(json::array({norm.ID(), self_coop_level, p.first}));

      /*
      double eq_c_level = EqCooperationLevel(norm);
      const double threshold = 0.2;
      if (eq_c_level > threshold) {
        // check local mutants as well
        auto p = MostRiskyMutant(norm, params);
        // std::pair<double,double> eq_c_level_local_eq0 = EqCooperationLevelWithLocalMutants(norm, params);
        output.push_back({norm.ID(), eq_c_level, p.first, p.second.ID()});
      }
      */
    }
    return output;
  };

  caravan::Option opt;
  opt.log_level = 2;
  opt.dump_log = "dump.log";
  caravan::Start(on_init, on_result_receive, do_task, MPI_COMM_WORLD, opt);

  if (my_rank == 0) {
    std::ofstream fout("results.txt");
    // sort results by cooperation level
    std::sort(results.begin(), results.end(), [](const auto& a, const auto& b) {
      return std::get<1>(a) + std::get<2>(a) > std::get<1>(b) + std::get<2>(b);
    });
    for (auto result: results) {
      fout << std::get<0>(result) << " " << std::get<1>(result) << " " << std::get<2>(result) << std::endl;
    }
    fout.close();
  }

  MPI_Finalize();

  return 0;
}
