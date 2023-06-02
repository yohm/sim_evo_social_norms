#include <iostream>
#include <thread>
#include <chrono>
#include <pagmo/algorithms/pso_gen.hpp>
#include <pagmo/problem.hpp>
#include <pagmo/batch_evaluators/member_bfe.hpp>
#include <mpi.h>
#include "Norm.hpp"
#include "PrivRepGame.hpp"


using namespace pagmo;

double EquilibriumCoopLevel(const std::array<double,20>& norm_vec) {
  EvolPrivRepGame::SimulationParameters params;
  params.n_init = 1e4;
  params.n_steps = 1e4;

  Norm norm = Norm::FromSerialized(norm_vec);

  EvolPrivRepGameAllCAllD evol(params, 5.0, 1.0);

  auto selfc_rho_eq = evol.EquilibriumCoopLevelAllCAllD(norm);
  double self_cooperation_level = std::get<0>(selfc_rho_eq);
  auto eq = std::get<2>(selfc_rho_eq);

  double c = self_cooperation_level * eq[0] + 1.0 * eq[1];
  IC(c);
  return c;
}

double EquilibriumCoopLevelWithoutR2(const std::array<double,12>& norm_vec) {
  EvolPrivRepGame::SimulationParameters params;
  params.n_init = 1e4;
  params.n_steps = 1e4;

  EvolPrivRepGameAllCAllD evol(params, 5.0, 1.0);

  AssessmentRule rd{{norm_vec[0], norm_vec[1], norm_vec[2], norm_vec[3], norm_vec[4], norm_vec[5], norm_vec[6], norm_vec[7]}};
  ActionRule ad{{norm_vec[8], norm_vec[9], norm_vec[10], norm_vec[11]}};
  Norm norm(rd, AssessmentRule::KeepRecipient(), ad);

  auto selfc_rho_eq = evol.EquilibriumCoopLevelAllCAllD(norm);
  double self_cooperation_level = std::get<0>(selfc_rho_eq);
  auto eq = std::get<2>(selfc_rho_eq);

  double c = self_cooperation_level * eq[0] + 1.0 * eq[1];
  IC(c);
  return c;
}

template <size_t input_dim>
struct my_problem {
  // Implementation of the objective function.
  vector_double fitness(const vector_double& dv) const;

  vector_double batch_fitness(const vector_double& dvs) const {
    size_t num_dvs = dvs.size() / input_dim;
    vector_double local_fitness(num_dvs, 0.0);

    int rank, num_procs;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // std::cerr << "batch_fitness called: " << rank << '/' << num_procs << std::endl;

    for (size_t i = rank; i < num_dvs; i+=num_procs) {
      size_t idx = i * input_dim;
      auto first = dvs.begin() + idx;
      auto last = dvs.begin() + idx + input_dim;
      std::vector<double> dv(first, last);
      double f = fitness(dv)[0];
      local_fitness[i] = f;
    }

    vector_double fitness(local_fitness.size(), 0.0);
    MPI_Allreduce(local_fitness.data(), fitness.data(), local_fitness.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    // std::cerr << "batch_fitness done: " << rank << '/' << num_procs << std::endl;

    return fitness;
  }

  // Implementation of the box bounds.
  std::pair<vector_double, vector_double> get_bounds() const
  {
    return {std::vector<double>(input_dim, 0.0), std::vector<double>(input_dim, 1.0)};
  }
};

template <>
vector_double my_problem<20>::fitness(const vector_double &dv) const {
  std::array<double,20> norm_vec{};
  for (int i = 0; i < 20; ++i) {
    norm_vec[i] = dv[i];
  }
  double coop_level = EquilibriumCoopLevel(norm_vec);
  return {-coop_level}; // pagmo solves minimization problem
}

template <>
vector_double my_problem<12>::fitness(const vector_double &dv) const {
  std::array<double,12> norm_vec{};
  for (int i = 0; i < 12; ++i) {
    norm_vec[i] = dv[i];
  }
  double coop_level = EquilibriumCoopLevelWithoutR2(norm_vec);
  return {-coop_level}; // pagmo solves minimization problem
}

int main(int argc, char* argv[]) {

  MPI_Init(&argc, &argv);
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  problem prob{my_problem<12>{}};
  // problem prob{my_problem<20>{}};

  // 2 - Instantiate a pagmo algorithm
  pagmo::random_device::set_seed(12345);
  unsigned max_gen = 100u;
  pso_gen pso_algo{max_gen};
  pso_algo.set_bfe(bfe{member_bfe{}});
  algorithm algo{pso_algo};
  if (my_rank == 0) algo.set_verbosity(max_gen / 10);

  // 3 - Instantiate a population of size 50
  population pop{prob, member_bfe{}, 50, 123456};

  // 4 - Let population evolve
  auto start = std::chrono::high_resolution_clock::now();
  pop = algo.evolve(pop);
  if( my_rank == 0) {
    std::cout << pop.champion_f()[0] << "\n";
    for (auto x : pop.champion_x()) {
      std::cout << x << ' ';
    }
  }
  auto end = std::chrono::high_resolution_clock::now();

  if (my_rank == 0) {
    std::cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms\n";
  }

  MPI_Finalize();

  return 0;
}

