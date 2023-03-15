#include <iostream>
#include <cassert>
#include <icecream.hpp>
#include <vector>
#include <random>
#include <mpi.h>
#include <caravan.hpp>
#include "Norm.hpp"
#include "PrivRepGame.hpp"


constexpr Reputation B = Reputation::B, G = Reputation::G;
constexpr Action C = Action::C, D = Action::D;

const EvolPrivRepGameAllCAllD* p_evol = nullptr;
const EvolPrivRepGameAllCAllD GetEvol() {
  if (p_evol) {
    return *p_evol;
  }
  else {
    EvolPrivRepGame::SimulationParameters params;
    params.n_init = 1e4;
    params.n_steps = 1e4;
    p_evol = new EvolPrivRepGameAllCAllD(30, params, 5.0, 1.0);
    return *p_evol;
  }
}

// return <ID, cooperation level> for a given norm
double EqCooperationLevel(const Norm& norm) {
  const auto evol = GetEvol();
  auto pc_rho_eq = evol.EquilibriumCoopLevelAllCAllD(norm);
  double self_cooperation_level = std::get<0>(pc_rho_eq);
  auto eq = std::get<2>(pc_rho_eq);
  return self_cooperation_level * eq[0] + 1.0 * eq[1];
}

void ComprehensiveSearchWithoutR2() {
  EvolPrivRepGameAllCAllD evol = GetEvol();

  std::vector< std::pair<int,double> > results;

  for (int j = 0; j < 256; j++) {
    std::cerr << "j = " << j << std::endl;
    AssessmentRule R1 = AssessmentRule::MakeDeterministicRule(j);
    AssessmentRule R2 = AssessmentRule::KeepRecipient();

    for (int i = 0; i < 16; i++) {
      ActionRule P = ActionRule::MakeDeterministicRule(i);
      Norm norm(R1, R2, P);
      if (norm.ID() < norm.SwapGB().ID()) {
        continue;
      }

      double c = EqCooperationLevel(norm);
      if (c > 0.3) {
        results.emplace_back(norm.ID(), c);
      }
    }
  }

  // sort results by cooperation level
  std::sort(results.begin(), results.end(), [](const auto& a, const auto& b) { return a.second > b.second; });
  // print top 10 from results
  for (auto result: results) {
    Norm norm = Norm::ConstructFromID(result.first);
    std::cout << norm.Inspect() << " " << result.second << std::endl;
  }
}


int main(int argc, char** argv) {
  // ComprehensiveSearchWithoutR2();

  MPI_Init(&argc, &argv);

  using namespace nlohmann;

  std::function<void(caravan::Queue&)> on_init = [](caravan::Queue& q) {
    const int i_max = 256, j_max = 256;
    // const int i_max = 60, j_max = 60;
    for (int i = 0; i < i_max; i++) {
      for (int j = 0; j < j_max; j++) {
        json input = {i, j};
        q.Push(input);
      }
    }
  };

  std::vector< std::pair<int,double> > results;
  std::function<void(int64_t, const json&, const json&, caravan::Queue&)> on_result_receive = [&results](int64_t task_id, const json& input, const json& output, caravan::Queue& q) {
    std::cerr << "task: " << task_id << " has finished: input: " << input << ", output: " << output << "\n";
    for (auto result: output) {
      results.emplace_back(result.at(0).get<int>(), result.at(1).get<double>());
    }
  };

  std::function<json(const json& input)> do_task = [](const json& input) {
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
      double eq_c_level = EqCooperationLevel(norm);
      double threshold = 0.2;
      if (eq_c_level > threshold) {
        output.push_back({norm.ID(), eq_c_level});
      }
    }
    return output;
  };

  caravan::Option opt;
  opt.log_level = 2;
  opt.dump_log = "dump.log";
  caravan::Start(on_init, on_result_receive, do_task, MPI_COMM_WORLD, opt);

  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  if (my_rank == 0) {
    std::ofstream fout("results.txt");
    // sort results by cooperation level
    std::sort(results.begin(), results.end(), [](const auto& a, const auto& b) { return a.second > b.second; });
    // print top 10 from results
    for (auto result: results) {
      Norm norm = Norm::ConstructFromID(result.first);
      fout << norm.Inspect() << " " << result.second << std::endl;
    }
    fout.close();
  }


  MPI_Finalize();

  return 0;
}
