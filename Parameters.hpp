//
// Created by Yohsuke Murase on 2023/05/23.
//

#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <nlohmann/json.hpp>
#include "PrivRepGame.hpp"

struct Parameters {
  size_t n_init;
  size_t n_steps;
  size_t N;
  double q;
  double mu_percept;
  double benefit;
  double beta;
  uint64_t seed;
  Parameters() : n_init(1e4), n_steps(1e4), N(30), q(0.9), mu_percept(0.05), benefit(5.0), beta(1.0), seed(123456789) {};

  EvolPrivRepGame::SimulationParameters ToEvolParams() const {
    return EvolPrivRepGame::SimulationParameters(N, n_init, n_steps, q, mu_percept, seed);
  }

  NLOHMANN_DEFINE_TYPE_INTRUSIVE_WITH_DEFAULT(Parameters, n_init, n_steps, N, q, mu_percept, benefit, beta, seed);
};

struct ParametersBatch {
  size_t n_init;
  size_t n_steps;
  size_t N;
  double q;
  double mu_percept;
  std::vector< std::pair<double,double> > benefit_beta_vec;
  uint64_t seed;
  ParametersBatch() : n_init(1e4), n_steps(1e4), N(30), q(0.9), mu_percept(0.05), seed(123456789) {
    benefit_beta_vec.emplace_back(5.0, 1.0);
  };

  Parameters ParameterAt(size_t i) const {  // return i-th parameter
    Parameters p;
    p.n_init = n_init;
    p.n_steps = n_steps;
    p.N = N;
    p.q = q;
    p.mu_percept = mu_percept;
    p.benefit = benefit_beta_vec[i].first;
    p.beta = benefit_beta_vec[i].second;
    p.seed = seed;
    return p;
  }

  EvolPrivRepGame::SimulationParameters ToEvolParams() const {
    return {N, n_init, n_steps, q, mu_percept, seed};
  }

  NLOHMANN_DEFINE_TYPE_INTRUSIVE_WITH_DEFAULT(ParametersBatch, n_init, n_steps, N, q, mu_percept, benefit_beta_vec, seed);
};

#endif //PARAMETERS_HPP
