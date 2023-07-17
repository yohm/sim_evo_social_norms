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
  double mu_assess;
  double benefit;
  double sigma_in_times_b;
  uint64_t seed;
  Parameters() : n_init(1e4), n_steps(1e4), N(30), q(1.0), mu_percept(0.0), mu_assess(0.05), benefit(5.0), sigma_in_times_b(5.0), seed(123456789) {};

  EvolPrivRepGame::SimulationParameters ToEvolParams() const {
    return EvolPrivRepGame::SimulationParameters(N, n_init, n_steps, q, mu_percept, mu_assess, seed);
  }
  double sigma_in() const {
    return sigma_in_times_b / (benefit-1.0);
  }

  NLOHMANN_DEFINE_TYPE_INTRUSIVE_WITH_DEFAULT(Parameters, n_init, n_steps, N, q, mu_percept, mu_assess, benefit, sigma_in_times_b, seed);
};

struct ParametersBatch {
  size_t n_init;
  size_t n_steps;
  size_t N;
  double q;
  double mu_percept;
  double mu_assess;
  private:
  std::vector< std::pair<double,double> > benefit_sigma_in_times_b_vec;
  public:
  std::vector<std::pair<double,double>> benefit_sigma_in_vec() const {
    std::vector<std::pair<double,double>> v;
    for (const auto& p : benefit_sigma_in_times_b_vec) {
      v.emplace_back(p.first, p.second / (p.first-1.0));
    }
    return v;
  }
  uint64_t seed;
  ParametersBatch() : n_init(1e4), n_steps(1e4), N(30), q(0.9), mu_percept(0.0), mu_assess(0.05), seed(123456789) {
    benefit_sigma_in_times_b_vec.emplace_back(5.0, 1.0);
  };

  Parameters ParameterAt(size_t i) const {  // return i-th parameter
    Parameters p;
    p.n_init = n_init;
    p.n_steps = n_steps;
    p.N = N;
    p.q = q;
    p.mu_percept = mu_percept;
    p.mu_assess = mu_assess;
    p.benefit = benefit_sigma_in_times_b_vec[i].first;
    p.sigma_in_times_b = benefit_sigma_in_times_b_vec[i].second;
    p.seed = seed;
    return p;
  }

  EvolPrivRepGame::SimulationParameters ToEvolParams() const {
    return {N, n_init, n_steps, q, mu_percept, mu_assess, seed};
  }

  NLOHMANN_DEFINE_TYPE_INTRUSIVE_WITH_DEFAULT(ParametersBatch, n_init, n_steps, N, q, mu_percept, mu_assess, benefit_sigma_in_times_b_vec, seed);
};

#endif //PARAMETERS_HPP
