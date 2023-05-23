#include <iostream>
#include <vector>
#include <random>
#include <Eigen/Dense>
#include <icecream.hpp>
#include "Norm.hpp"

class PrivateRepGame {
public:
  using population_t = std::vector<std::pair<Norm, size_t>>;  // vector of StrategyID & its size

  PrivateRepGame(const population_t& pop, uint64_t _seed) : population(pop), rnd(_seed), uni(0.0, 1.0) {
    for (size_t i = 0; i < population.size(); i++) {
      auto kv = population[i];
      for (size_t ii = 0; ii < kv.second; ii++) {
        norms.emplace_back(kv.first);
        norm_index.emplace_back(i);
      }
    }
    N = norms.size();
    M.assign(N*N, Reputation::G);
    // initialize coop_count with NxN matrix of 0
    coop_count.assign(N*N, 0);
    game_count.assign(N*N, 0);
    good_count.assign(N*N, 0);
  }

  // t_max : number of steps
  // q : observation probability
  // epsilon :
  void Update(size_t t_max, double q, double mu_percept, bool count_good = true) {
    for (size_t t = 0; t < t_max; t++) {
      // randomly choose donor & recipient
      size_t donor = static_cast<size_t>(R01() * N);
      size_t recip = (donor + static_cast<size_t>(R01() * (N-1)) + 1) % N;
      assert(donor != recip);

      double c_prob = norms[donor].P.CProb(M[donor*N+donor], M[donor*N+recip]);
      Action A = (c_prob == 1.0 || R01() < c_prob) ? Action::C : Action::D;

      if (A == Action::C) {
        coop_count[donor*N+recip]++;
      }
      game_count[donor*N+recip]++;

      // updating the images from observers' viewpoint
      for (size_t obs = 0; obs < N; obs++) {
        if (obs == donor || obs == recip || R01() < q) {  // observe with probability q
          Action a_obs = A;
          if (mu_percept > 0.0 && R01() < mu_percept) {
            a_obs = FlipAction(A);
          }

          // update donor's reputation
          double g_prob_donor = norms[obs].Rd.GProb(M[obs*N+donor], M[obs*N+recip], a_obs);
          if (g_prob_donor == 1.0) {
            M[obs*N+donor] = Reputation::G;
          }
          else if (g_prob_donor == 0.0) {
            M[obs*N+donor] = Reputation::B;
          }
          else {
            M[obs*N+donor] = (R01() < g_prob_donor) ? Reputation::G : Reputation::B;
          }

          // update recipient's reputation
          double g_prob_recip = norms[obs].Rr.GProb(M[obs*N+donor], M[obs*N+recip], a_obs);
          if (g_prob_recip == 1.0) {
            M[obs*N+recip] = Reputation::G;
          }
          else if (g_prob_recip == 0.0) {
            M[obs*N+recip] = Reputation::B;
          }
          else {
            M[obs*N+recip] = (R01() < g_prob_recip) ? Reputation::G : Reputation::B;
          }
        }
      }

      // count good
      if (count_good) {
        for (size_t i = 0; i < N; i++) {
          for (size_t j = 0; j < N; j++) {
            if (M[i*N+j] == Reputation::G) {
              good_count[i*N+j]++;
            }
          }
        }
        total_update++;
      }

    }
  }

  using count_t = std::vector<size_t>;

  // number of C between i-donor and j-recipient
  count_t CoopCount() const { return coop_count; }

  // number of games between i-donor and j-recipient
  count_t GameCount() const { return game_count; }

  // system-wide cooperation level
  double SystemWideCooperationLevel() const {
    double coop = 0.0;
    double total = 0.0;
    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < N; j++) {
        coop += coop_count[i*N+j];
        total += game_count[i*N+j];
      }
    }
    return coop / total;
  }

  // donation levels of individuals
  // <probability of receiving benefit, probability of paying cost>
  // payoff
  std::vector<std::pair<double,double>> IndividualCooperationLevels() const {
    std::vector<std::pair<double,double>> coop_rates(N, {0.0, 0.0});
    for (size_t i = 0; i < N; i++) {
      double coop_total_out = 0.0, game_total_out = 0.0;
      for (size_t j = 0; j < N; j++) {
        if (i == j) {
          continue;
        }
        coop_total_out += coop_count[i*N+j];
        game_total_out += game_count[i*N+j];
      }
      double coop_total_in = 0.0, game_total_in = 0.0;
      for (size_t j = 0; j < N; j++) {
        if (i == j) {
          continue;
        }
        coop_total_in += coop_count[j*N+i];
        game_total_in += game_count[j*N+i];
      }
      coop_rates[i] = std::make_pair(coop_total_in / game_total_in, coop_total_out / game_total_out);
    }
    return coop_rates;
  }

  // norm-wise cooperation levels
  // return: vector c_levels
  //   c_levels[i][j] : cooperation level of i-th norm toward j-th norm
  std::vector<std::vector<double>> NormCooperationLevels() const {
    size_t n_norms = population.size();
    std::vector<std::vector<size_t>> coop_by_norm(n_norms, std::vector<size_t>(n_norms, 0) );
    std::vector<std::vector<size_t>> total_by_norm(n_norms, std::vector<size_t>(n_norms, 0) );
    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < N; j++) {
        size_t i_norm = norm_index[i];
        size_t j_norm = norm_index[j];
        coop_by_norm[i_norm][j_norm] += coop_count[i*N+j];
        total_by_norm[i_norm][j_norm] += game_count[i*N+j];
      }
    }
    std::vector<std::vector<double>> c_levels(n_norms, std::vector<double>(n_norms, 0.0));
    for (size_t i_norm = 0; i_norm < n_norms; i_norm++) {
      for (size_t j_norm = 0; j_norm < n_norms; j_norm++) {
        c_levels[i_norm][j_norm] = static_cast<double>(coop_by_norm[i_norm][j_norm]) / static_cast<double>(total_by_norm[i_norm][j_norm]);
      }
    }
    return c_levels;
  }

  // norm-wise good count
  // return: vector average reputation
  //   c_levels[i][j] : average reputation of j-th norm from the viewpoint of i-th norm
  std::vector<std::vector<double>> NormAverageReputation() const {
    size_t n_norms = population.size();
    std::vector<std::vector<double>> avg_rep(n_norms, std::vector<double>(n_norms, 0.0));

    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < N; j++) {
        size_t i_norm = norm_index[i];
        size_t j_norm = norm_index[j];
        avg_rep[i_norm][j_norm] += good_count[i*N+j];
      }
    }

    for (size_t i_norm = 0; i_norm < n_norms; i_norm++) {
      for (size_t j_norm = 0; j_norm < n_norms; j_norm++) {
        size_t c = population[i_norm].second * population[j_norm].second * total_update;
        avg_rep[i_norm][j_norm] /= static_cast<double>(c);
      }
    }

    return avg_rep;
  }

  void ResetCounts() {
    // reset coop_count and fill with 0
    std::fill(coop_count.begin(), coop_count.end(), 0);
    std::fill(game_count.begin(), game_count.end(), 0);
    std::fill(good_count.begin(), good_count.end(), 0);
    total_update = 0;
  }

  // print image matrix
  void PrintImage(std::ostream& out) const {
    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < N; j++) {
        out << ((M[i*N+j]==Reputation::G)?'.':'x');
      }
      out << "\n";
    }
  }

  population_t Population() const {
    return population;
  }

private:
  const population_t population;
  size_t N;
  std::mt19937_64 rnd;
  std::uniform_real_distribution<double> uni;
  std::vector<Norm> norms;
  std::vector<size_t> norm_index;  // [0,0,...,0, 1,1,....,1, 2,2,......,2]
  std::vector<Reputation> M;
  count_t coop_count;  // number of C between i-donor and j-recipient
  count_t game_count;  // number of games between i-donor and j-recipient
  count_t good_count;  // number of good reputations
  size_t total_update; // number of time steps in total
  double R01() { return uni(rnd); }
};


class EvolPrivRepGame {
public:
  struct SimulationParameters {
    SimulationParameters(size_t N = 50, size_t n_init = 1e6, size_t n_steps = 1e6, double q = 0.9, double mu_percept = 0.05, uint64_t seed = 123456789ull) :
        N(N), n_init(n_init), n_steps(n_steps), q(q), mu_percept(mu_percept), seed(seed) {};
    size_t N;  // population size
    size_t n_init, n_steps;  // simulation duration
    double q;  // observation probability
    double mu_percept;  // perception error probability
    uint64_t seed;  // random number seed
  };

  using norms_t = std::vector<Norm>;

  EvolPrivRepGame(const SimulationParameters& sim_param) : param(sim_param) {};
  SimulationParameters param;

  // rho[i][j] = fixation probability of a j-mutant into i-resident community
  std::vector<std::vector<double>> FixationProbabilities(const norms_t& norms, double benefit, double beta) const {
    size_t num_norms = norms.size();
    std::vector<std::vector<double>> rho(num_norms, std::vector<double>(num_norms, 0.0));

    for (size_t i = 0; i < num_norms; i++) {
      for (size_t j = i+1; j < num_norms; j++) {
        auto rho_ij_ji = FixationProbability(norms[i], norms[j], benefit, beta);
        rho[i][j] = rho_ij_ji.first;
        rho[j][i] = rho_ij_ji.second;
      }
    }

    return rho;
  }

  // self cooperation level
  double SelfCooperationLevel(const Norm& norm) const {
    PrivateRepGame game({{norm, param.N}}, param.seed);
    game.Update(param.n_init, param.q, param.mu_percept, false);
    game.ResetCounts();
    game.Update(param.n_steps, param.q, param.mu_percept, false);
    double self_coop_level = game.NormCooperationLevels()[0][0];
    return self_coop_level;
  }

  // first: fixation probability of resident j against resident i
  //        i.e., the probability to change from i to j
  // second: fixation probability of resident i against resident j
  //        i.e., the probability to change from j to i
  std::pair<double,double> FixationProbability(const Norm& norm_i, const Norm& norm_j, double benefit, double beta) const {
    const size_t N = param.N;
    std::vector<double> pi_i(N);  // pi_i[l]: payoff of resident i when l mutants exist
    std::vector<double> pi_j(N);  // pi_j[l]: payoff of mutant j when l mutants exist

    #pragma omp parallel for schedule(dynamic) default(none), shared(pi_i, pi_j, norm_i, norm_j, benefit, beta, N)
    for (size_t l = 1; l < N; l++) {
      PrivateRepGame game({{norm_i, N-l}, {norm_j, l}}, param.seed + l);
      game.Update(param.n_init, param.q, param.mu_percept, false);
      game.ResetCounts();
      game.Update(param.n_steps, param.q, param.mu_percept, false);
      auto coop_levles = game.IndividualCooperationLevels();
      double payoff_i_total = 0.0;
      for (size_t k = 0; k < N-l; k++) {
        payoff_i_total += benefit * coop_levles[k].first - coop_levles[k].second;
      }
      pi_i[l] = payoff_i_total / (N-l);
      double payoff_j_total = 0.0;
      for (size_t k = N-l; k < N; k++) {
        payoff_j_total += benefit * coop_levles[k].first - coop_levles[k].second;
      }
      pi_j[l] = payoff_j_total / l;
      // IC(l, pi_i[l], pi_j[l]);
    }

    double rho_1_inv = 1.0;
    // p_ij = 1 / (1 + sum_{l' != 1}^{N-1}  prod_{l=1}^{l_prime} exp{-beta * (pi_j[l] - pi_i[l]) }
    for (size_t l_prime = 1; l_prime < N; l_prime++) {
      double prod = 1.0;
      for (size_t l = 1; l <= l_prime; l++) {
        prod *= exp(-beta * (pi_j[l] - pi_i[l]));
      }
      rho_1_inv += prod;
    }
    double rho_1 = 1.0 / rho_1_inv;

    double rho_2_inv = 1.0;
    for (size_t l_prime = 1; l_prime < N; l_prime++) {
      double prod = 1.0;
      for (size_t l = 1; l <= l_prime; l++) {
        prod *= exp(-beta * (pi_i[N-l] - pi_j[N-l]));
      }
      rho_2_inv += prod;
    }
    double rho_2 = 1.0 / rho_2_inv;
    return std::make_pair(rho_1, rho_2);
  }

  // first: fixation probability of resident j against resident i
  //        i.e., the probability to change from i to j
  // second: fixation probability of resident i against resident j
  //        i.e., the probability to change from j to i
  std::vector<std::pair<double,double>> FixationProbabilityBatch(const Norm& norm_i, const Norm& norm_j, const std::vector<std::pair<double,double>>& benefit_beta_pairs) const {
    std::vector<std::pair<double,double>> ans;
    for (auto& [benefit,beta] : benefit_beta_pairs) {
      ans.push_back( FixationProbability(norm_i, norm_j, benefit, beta) );
    }
    return ans;
  }

  // arg: transition probability matrix p
  //      p[i][j] : fixation probability of a j-mutant into i-resident
  // return: stationary distribution vector
  static std::vector<double> EquilibriumPopulationLowMut(const std::vector<std::vector<double>>& fixation_probs) {
    size_t num_norms = fixation_probs.size();

    // construct transition matrix
    Eigen::MatrixXd T(num_norms, num_norms);

    for (size_t i = 0; i < num_norms; i++) {
      double sum = 0.0;
      for (size_t j = 0; j < num_norms; j++) {
        if (i == j) continue;
        T(j,i) = fixation_probs[i][j] / (1.0 - num_norms);
        sum += T(j, i);
      }
      T(i,i) = - sum;  // T(i,i) = 1.0 - sum; but we subtract I to calculate stationary distribution
    }
    for (size_t i = 0; i < num_norms; i++) {  // normalization condition
      T(num_norms-1, i) += 1.0;
    }

    Eigen::VectorXd b(num_norms);
    for (int i = 0; i < num_norms-1; i++) { b(i) = 0.0; }
    b(num_norms-1) = 1.0;

    Eigen::VectorXd x = T.colPivHouseholderQr().solve(b);

    std::vector<double> ans(num_norms, 0.0);
    for (int i = 0; i < num_norms; i++) {
      ans[i] = x(i);
    }
    return ans;
  }

};


// Evolutionary game between X and AllC and AllD
class EvolPrivRepGameAllCAllD {
public:
  EvolPrivRepGameAllCAllD(const EvolPrivRepGame::SimulationParameters& sim_param, double benefit, double selection_strength) :
    benefit(benefit), selection_strength(selection_strength), param(sim_param) {
    auto p = FixationProbsBetweenAllCAllD();
    rho_allc = p.first;
    rho_alld = p.second;
  }

  const double benefit;
  const double selection_strength;
  const EvolPrivRepGame::SimulationParameters param;
  double rho_allc, rho_alld;

  // fixation probability between AllC and AllD
  std::pair<double,double> FixationProbsBetweenAllCAllD() const {
    const size_t N = param.N;
    const Norm allc = Norm::AllC(), alld = Norm::AllD();

    // payoff when there are i AllC and N-i AllD players
    std::vector<double> pi_allc(N, 0.0), pi_alld(N, 0.0);

    #pragma omp parallel for schedule(dynamic) default(none), shared(pi_allc, pi_alld, allc, alld, N)
    for (size_t i = 1; i < N; i++) {
      // i AllC vs N-i AllD
      PrivateRepGame game({{allc, i}, {alld, N-i}}, param.seed);
      game.Update(param.n_init, param.q, param.mu_percept, false);
      game.ResetCounts();
      game.Update(param.n_steps, param.q, param.mu_percept, false);
      auto coop_levels = game.IndividualCooperationLevels();
      double payoff_i_total = 0.0;
      for (size_t k = 0; k < i; k++) {
        payoff_i_total += benefit * coop_levels[k].first - coop_levels[k].second;
      }
      pi_allc[i] = payoff_i_total / i;
      double payoff_j_total = 0.0;
      for (size_t k = i; k < N; k++) {
        payoff_j_total += benefit * coop_levels[k].first - coop_levels[k].second;
      }
      pi_alld[i] = payoff_j_total / (N-i);
    }
    // IC(pi_allc, pi_alld);

    // the probability that X is fixed in a population of AllC
    // p_ij = 1 / (1 + sum_{l' != 1}^{N-1}  prod_{l=1}^{l_prime} exp{-beta * (pi_j[l] - pi_i[l]) }
    double rho_inv = 1.0;
    for (size_t l_prime = 1; l_prime < N; l_prime++) {
      double prod = 1.0;
      for (size_t l = 1; l <= l_prime; l++) {
        prod *= exp(-selection_strength * (pi_allc[l] - pi_alld[l]));
      }
      rho_inv += prod;
    }
    double rho_allc = 1.0 / rho_inv;

    rho_inv = 1.0;
    for (size_t l_prime = 1; l_prime < N; l_prime++) {
      double prod = 1.0;
      for (size_t l = 1; l <= l_prime; l++) {
        prod *= exp(-selection_strength * (pi_alld[N-l] - pi_allc[N-l]));
      }
      rho_inv += prod;
    }
    double rho_alld = 1.0 / rho_inv;

    return std::make_pair(rho_allc, rho_alld);
  }

  // return self-cooperation-level, rho, equilibrium population
  std::tuple<double,std::vector<std::vector<double>>,std::vector<double>> EquilibriumCoopLevelAllCAllD(const Norm& norm) const {
    const size_t N = param.N;
    std::vector<double> pi_x_allc(N, 0.0), pi_allc(N, 0.0), pi_x_alld(N, 0.0), pi_alld(N, 0.0);
    double self_coop_level;
    // pi_x_allc[l] : payoff of X against AllC when there are l AllC and (N-l) residents
    // pi_allc[l]   : payoff of AllC against X when there are l AllC and (N-l) residents
    // pi_x_alld[l] : payoff of X against AllD when there are l AllD and (N-l) residents
    // pi_alld[l]   : payoff of AllD against X when there are l AllD and (N-l) residents

    #pragma omp parallel for schedule(dynamic) default(none), shared(pi_x_allc, pi_allc, pi_x_alld, pi_alld, norm, self_coop_level, N)
    for (size_t i = 0; i < 2*N-1; i++) {
      if (i == 0) {  // monomorphic population of X
        PrivateRepGame game({{norm, N}}, param.seed);
        game.Update(param.n_init, param.q, param.mu_percept, false);
        game.ResetCounts();
        game.Update(param.n_steps, param.q, param.mu_percept, false);
        self_coop_level = game.SystemWideCooperationLevel();
        pi_x_allc[0] = pi_x_alld[0] = self_coop_level * (benefit - 1.0);
      }
      else if (i < N) {
        size_t l = i;
        // (N-l) X vs l AllC
        PrivateRepGame game({{norm, N-l}, {Norm::AllC(), l}}, param.seed);
        game.Update(param.n_init, param.q, param.mu_percept, false);
        game.ResetCounts();
        game.Update(param.n_steps, param.q, param.mu_percept, false);
        auto coop_levels = game.IndividualCooperationLevels();
        double payoff_i_total = 0.0;
        for (size_t k = 0; k < N-l; k++) {
          payoff_i_total += benefit * coop_levels[k].first - coop_levels[k].second;
        }
        pi_x_allc[l] = payoff_i_total / (N-l);
        double payoff_j_total = 0.0;
        for (size_t k = N-l; k < N; k++) {
          payoff_j_total += benefit * coop_levels[k].first - coop_levels[k].second;
        }
        pi_allc[l] = payoff_j_total / l;
      }
      else if (i >= N) {
        // (N-l) X vs l AllD
        size_t l = i - N + 1;
        PrivateRepGame game({{norm, N-l}, {Norm::AllD(), l}}, param.seed);
        game.Update(param.n_init, param.q, param.mu_percept, false);
        game.ResetCounts();
        game.Update(param.n_steps, param.q, param.mu_percept, false);
        auto coop_levels = game.IndividualCooperationLevels();
        double payoff_i_total = 0.0;
        for (size_t k = 0; k < N-l; k++) {
          payoff_i_total += benefit * coop_levels[k].first - coop_levels[k].second;
        }
        pi_x_alld[l] = payoff_i_total / (N-l);
        double payoff_j_total = 0.0;
        for (size_t k = N-l; k < N; k++) {
          payoff_j_total += benefit * coop_levels[k].first - coop_levels[k].second;
        }
        pi_alld[l] = payoff_j_total / l;
      }
    }
    // IC(pi_x_allc, pi_allc, pi_x_alld, pi_alld);

    // calculate the fixation probabilities
    std::vector<std::vector<double>> rho(3, std::vector<double>(3, 0.0));

    {
      // the probability that X is fixed in a population of AllC
      // p_ij = 1 / (1 + sum_{l' != 1}^{N-1}  prod_{l=1}^{l_prime} exp{-beta * (pi_j[l] - pi_i[l]) }
      double rho_inv = 1.0;
      for (size_t l_prime = 1; l_prime < N; l_prime++) {
        double prod = 1.0;
        for (size_t l = 1; l <= l_prime; l++) {
          prod *= exp(-selection_strength * (pi_allc[l] - pi_x_allc[l]));
        }
        rho_inv += prod;
      }
      rho[0][1] = 1.0 / rho_inv;
    }
    {
      // the probability that X is fixed in a population of AllD
      double rho_inv = 1.0;
      for (size_t l_prime = 1; l_prime < N; l_prime++) {
        double prod = 1.0;
        for (size_t l = 1; l <= l_prime; l++) {
          prod *= exp(-selection_strength * (pi_alld[l] - pi_x_alld[l]));
        }
        rho_inv += prod;
      }
      rho[0][2] = 1.0 / rho_inv;
    }
    {
      // the probability that AllC is fixed in a population of X
      double rho_inv = 1.0;
      for (size_t l_prime = 1; l_prime < N; l_prime++) {
        double prod = 1.0;
        for (size_t l = 1; l <= l_prime; l++) {
          prod *= exp(-selection_strength * (pi_x_allc[N-l] - pi_allc[N-l]));
        }
        rho_inv += prod;
      }
      rho[1][0] = 1.0 / rho_inv;
    }
    {
      // the probability that AllD is fixed in a population of X
      double rho_inv = 1.0;
      for (size_t l_prime = 1; l_prime < N; l_prime++) {
        double prod = 1.0;
        for (size_t l = 1; l <= l_prime; l++) {
          prod *= exp(-selection_strength * (pi_x_alld[N-l] - pi_alld[N-l]));
        }
        rho_inv += prod;
      }
      rho[2][0] = 1.0 / rho_inv;
    }
    rho[1][2] = rho_alld;
    rho[2][1] = rho_allc;

    // calculate the equilibrium populations
    auto eq = EvolPrivRepGame::EquilibriumPopulationLowMut(rho);

    return std::make_tuple(self_coop_level, rho, eq);
  };
};


class EvolPrivRepGameFiniteMutationRateAllCAllD {
public:
  using SimulationParameters = EvolPrivRepGame::SimulationParameters;

  EvolPrivRepGameFiniteMutationRateAllCAllD(const Norm &norm, const SimulationParameters &sim_param) :
      norm(norm), param(sim_param) {
    games.resize(param.N+1);
    #pragma omp parallel for schedule(dynamic,1) shared(games)
    for (size_t nf = 0; nf <= param.N; nf++) {
      if (nf % 5 == 0) { std::cerr << "simulating games at nf=" << nf << std::endl; }
      for (size_t nc = 0; nc <= param.N-nf; nc++) {
        games[nf].push_back(RunSimulationAt(nf, nc));
      }
    }
  };

  const Norm norm;
  SimulationParameters param;
  std::vector<std::vector<PrivateRepGame>> games;

  struct EquilibriumState {
    std::vector<std::vector<double>> frequency;
    std::vector<std::vector<double>> cooperation_level;  // cooperation level at state[nf][nc]
    const size_t N;
    EquilibriumState(size_t _N) : N(_N) {
      frequency.resize(N+1);
      cooperation_level.resize(N+1);
      for (size_t i = 0; i <= N; i++) {
        frequency[i].resize(N+1-i, 0.0);
        cooperation_level[i].resize(N+1-i, 0.0);
      }
    }
    double OverallCooperationLevel() const {
      double total = 0.0;
      for (size_t nf = 0; nf < frequency.size(); nf++) {
        for (size_t nc = 0; nc < frequency[nf].size(); nc++) {
          total += frequency[nf][nc] * cooperation_level[nf][nc];
        }
      }
      return total;
    }
    std::array<double,3> OverallAbundances() const {
      std::array<double,3> abundances = {0.0, 0.0, 0.0};
      for (size_t nf = 0; nf < frequency.size(); nf++) {
        for (size_t nc = 0; nc < frequency[nf].size(); nc++) {
          size_t nd = N - nf - nc;
          abundances[0] += frequency[nf][nc] * nf / (double)N;
          abundances[1] += frequency[nf][nc] * nc / (double)N;
          abundances[2] += frequency[nf][nc] * nd / (double)N;
        }
      }
      return abundances;
    }
  };

  EquilibriumState CalculateEquilibrium(double benefit, double beta, double mu) const {
    // [TODO] ideally, we should calculate the average of transition probabilities, not payoffs

    // index of (nf, nc) in the one-dimensional vector
    std::vector<std::vector<size_t>> pair2index;
    std::vector<std::pair<size_t,size_t>> index2pair;
    for (size_t nf = 0; nf <= param.N; nf++) {
      pair2index.emplace_back();
      for (size_t nc = 0; nc <= param.N-nf; nc++) {
        index2pair.emplace_back(nf, nc);
        pair2index[nf].emplace_back(index2pair.size()-1);
      }
    }
    // IC(pair2index, index2pair);

    EquilibriumState equilibrium_state(param.N);

    size_t num_states = index2pair.size();
    Eigen::MatrixXd W(num_states, num_states);
    W.setZero();

    // calculate the transition probabilities
    for (size_t nf = 0; nf <= param.N; nf++) {
      for (size_t nc = 0; nc <= param.N - nf; nc++) {
        size_t nd = param.N - nf - nc;
        const PrivateRepGame g = games.at(nf).at(nc);
        equilibrium_state.cooperation_level[nf][nc] = g.SystemWideCooperationLevel();
        auto [w_f_c, w_f_d, w_c_f, w_c_d, w_d_f, w_d_c] = CalculateTransitionProbabilitiesFromGame(benefit, beta, mu, g);
        if (nf > 0) {
          W(pair2index[nf-1][nc+1], pair2index[nf][nc]) = w_f_c;
          W(pair2index[nf-1][nc]  , pair2index[nf][nc]) = w_f_d;
        }
        if (nc > 0) {
          W(pair2index[nf+1][nc-1], pair2index[nf][nc]) = w_c_f;
          W(pair2index[nf][nc-1]  , pair2index[nf][nc]) = w_c_d;
        }
        if (nd > 0) {
          W(pair2index[nf+1][nc]  , pair2index[nf][nc]) = w_d_f;
          W(pair2index[nf][nc+1]  , pair2index[nf][nc]) = w_d_c;
        }
        W(pair2index[nf][nc], pair2index[nf][nc]) = - w_f_c - w_f_d - w_c_f - w_c_d - w_d_f - w_d_c;
        // W(i,i) = 1.0 - sum; but we subtract I to calculate stationary distribution
      }
    }
    for (size_t i = 0; i < num_states; i++) {
      W(num_states-1, i) += 1.0;
      // normalization condition of the answer
    }

    Eigen::VectorXd b(num_states);
    b.setZero();
    b(num_states-1) = 1.0;

    // solve (W-I)x = 0
    Eigen::VectorXd x = W.colPivHouseholderQr().solve(b);

    std::vector<double> ans(num_states, 0.0);
    for (int i = 0; i < num_states; i++) {
      ans[i] = x(i);
    }

    // calculate the equilibrium state
    for (size_t i = 0; i < num_states; i++) {
      auto [nf, nc] = index2pair[i];
      equilibrium_state.frequency[nf][nc] = ans[i];
    }

    return equilibrium_state;
  }

private:
  // run simulation at [nf, nc, N-nf-nc] and return the game
  PrivateRepGame RunSimulationAt(size_t nf, size_t nc) const {
    PrivateRepGame game({{norm, nf}, {Norm::AllC(), nc}, {Norm::AllD(), param.N-nf-nc}}, param.seed);
    game.Update(param.n_init, param.q, param.mu_percept, false);
    game.ResetCounts();
    game.Update(param.n_steps, param.q, param.mu_percept, true);
    return game;
  }

  // calculate payoffs <pi_f, pi_c, pi_d>
  std::array<double,3> CalculatePayoffsFromGame(double benefit, const PrivateRepGame& game) const {
    auto coop_levels = game.IndividualCooperationLevels();
    auto pop = game.Population();
    size_t nf = pop.at(0).second;  // number of the focal norm players
    size_t nc = pop.at(1).second;  // number of AllC players
    size_t nd = pop.at(2).second;  // number of AllD players

    double pi_f = std::numeric_limits<double>::signaling_NaN();
    double pi_c = std::numeric_limits<double>::signaling_NaN();
    double pi_d = std::numeric_limits<double>::signaling_NaN();
    if (nf > 0) {
      pi_f = 0.0;
      for (size_t i = 0; i < nf; i++) {
        pi_f += benefit * coop_levels[i].first - coop_levels[i].second;
      }
      pi_f /= nf;
    }
    if (nc > 0) {
      pi_c = 0.0;
      for (size_t i = nf; i < nf+nc; i++) {
        pi_c += benefit * coop_levels[i].first - coop_levels[i].second;
      }
      pi_c /= nc;
    }
    if (nd > 0) {
      pi_d = 0.0;
      for (size_t i = nf+nc; i < param.N; i++) {
        pi_d += benefit * coop_levels[i].first - coop_levels[i].second;
      }
      pi_d /= nd;
    }
    return {pi_f, pi_c, pi_d};
  }

  // calculate the transition probabilities at [nf, nc, N-nf-nc]
  // return: [w(f->c), w(f->d), w(c->f), w(c->d), w(d->f), w(d->c)]
  std::array<double,6> CalculateTransitionProbabilitiesFromGame(double benefit, double beta, double mu, const PrivateRepGame& game) const {
    auto [pi_f, pi_c, pi_d] = CalculatePayoffsFromGame(benefit, game);

    auto pop = game.Population();
    size_t nf = pop.at(0).second;  // number of the focal norm players
    size_t nc = pop.at(1).second;  // number of AllC players
    size_t nd = pop.at(2).second;  // number of AllD players
    const size_t N = param.N;

    double rho_f = nf / (double) N;
    double rho_c = nc / (double) N;
    double rho_d = nd / (double) N;

    // transition probability from f to c
    // rho_f * { mu/(2) + (1-mu) * rho_c * [1 + exp(-beta * (pi_c - pi_f)) ]^-1 }
    double w_f_c = 0.0;
    if (nf > 0) {
      w_f_c = rho_f * mu * 0.5;  // flow f->c via mutation
      if (nc > 0) {  // flow f->c via imitation
        double P_f_c = 1.0 / (1.0 + exp(-beta * (pi_c - pi_f)));
        w_f_c += rho_f * (1.0 - mu) * (nc/(double)(N-1)) * P_f_c;
      }
    }

    // transition probability from f to d
    double w_f_d = 0.0;
    if (nf > 0) {
      w_f_d = rho_f * mu * 0.5;  // flow f->d via mutation
      if (nd > 0) {  // flow f->d via imitation
        double P_f_d = 1.0 / (1.0 + exp(-beta * (pi_d - pi_f)));
        w_f_d += rho_f * (1.0 - mu) * (nd/(double)(N-1)) * P_f_d;
      }
    }

    // transition probability from c to f
    double w_c_f = 0.0;
    if (nc > 0) {
      w_c_f = rho_c * mu * 0.5;  // flow c->f via mutation
      if (nf > 0) {  // flow c->f via imitation
        double P_c_f = 1.0 / (1.0 + exp(-beta * (pi_f - pi_c)));
        w_c_f += rho_c * (1.0 - mu) * (nf/(double)(N-1)) * P_c_f;
      }
    }

    // transition probability from c to d
    double w_c_d = 0.0;
    if (nc > 0) {
      w_c_d = rho_c * mu * 0.5;  // flow c->d via mutation
      if (nd > 0) {  // flow c->d via imitation
        double P_c_d = 1.0 / (1.0 + exp(-beta * (pi_d - pi_c)));
        w_c_d += rho_c * (1.0 - mu) * (nd/(double)(N-1)) * P_c_d;
      }
    }

    // transition probability from d to f
    double w_d_f = 0.0;
    if (nd > 0) {
      w_d_f = rho_d * mu * 0.5;  // flow d->f via mutation
      if (nf > 0) {  // flow d->f via imitation
        double P_d_f = 1.0 / (1.0 + exp(-beta * (pi_f - pi_d)));
        w_d_f += rho_d * (1.0 - mu) * (nf/(double)(N-1)) * P_d_f;
      }
    }

    // transition probability from d to c
    double w_d_c = 0.0;
    if (nd > 0) {
      w_d_c = rho_d * mu * 0.5;  // flow d->c via mutation
      if (nc > 0) {  // flow d->c via imitation
        double P_d_c = 1.0 / (1.0 + exp(-beta * (pi_c - pi_d)));
        w_d_c += rho_d * (1.0 - mu) * (nc/(double)(N-1)) * P_d_c;
      }
    }

    return {w_f_c, w_f_d, w_c_f, w_c_d, w_d_f, w_d_c};
  }
};