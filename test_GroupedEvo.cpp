#include <iostream>
#include <fstream>
#include <chrono>
#include <regex>
#include <icecream.hpp>
#include "GroupedEvo.hpp"

#include <gtest/gtest.h>


TEST(GroupedEvo, InterGroupImitationProb) {
  double p_imit = GroupedEvo::InterGroupImitationProb(1.0, 0.5, 5.0, 1.0);
  double expected = 1.0 / (1.0 + std::exp(1.0 * (4.0 - 2.0)));
  EXPECT_NEAR(p_imit, expected, 0.002);

  double p_imit2 = GroupedEvo::InterGroupImitationProb(1.0, 0.25, 5.0, 0.7);
  double expected2 = 1.0 / (1.0 + std::exp(0.7 * (4.0 - 1.0)));
  EXPECT_NEAR(p_imit2, expected2, 0.002);
}

TEST(GroupedEvo, CalcAlpha) {
  std::vector<size_t> norm_ids = {0, 1, 2};
  std::vector<std::vector<double>> p_fix = {{1.0/3.0, 0.5, 1.0}, {0.25, 1.0/3.0, 0.1}, {0.0, 0.5, 0.0}};
  std::vector<double> self_coop_levels = {0.0, 0.5, 1.0};
  GroupedEvo ge(norm_ids, p_fix, self_coop_levels);
  // 0-1
  {
    double alpha = ge.CalcAlpha(0, 1, 5.0, 1.0);
    double plus = 0.5 / (1.0 + std::exp(1.0 * (0.0 - 2.0)));
    double minus = 0.25 / (1.0 + std::exp(1.0 * (2.0 - 0.0)));
    double expected = plus - minus;
    EXPECT_NEAR(alpha, expected, 0.002);

    // alpha should be anti-symmetric
    double alpha_inv = ge.CalcAlpha(1, 0, 5.0, 1.0);
    EXPECT_NEAR(alpha_inv, -expected, 0.002);  //
  }
  // 0-2
  {
    double alpha = ge.CalcAlpha(0, 2, 5.0, 1.0);
    double plus = 1.0 / (1.0 + std::exp(1.0 * (0.0 - 4.0)));
    double minus = 0.0 / (1.0 + std::exp(1.0 * (4.0 - 0.0)));
    double expected = plus - minus;
    EXPECT_NEAR(alpha, expected, 0.002);

    // alpha should be anti-symmetric
    double alpha_inv = ge.CalcAlpha(2, 0, 5.0, 1.0);
    EXPECT_NEAR(alpha_inv, -expected, 0.002);  //
  }
  // 1-2
  {
    double alpha = ge.CalcAlpha(1, 2, 5.0, 1.0);
    double plus = 0.1 / (1.0 + std::exp(1.0 * (2.0 - 4.0)));
    double minus = 0.5 / (1.0 + std::exp(1.0 * (4.0 - 2.0)));
    double expected = plus - minus;
    EXPECT_NEAR(alpha, expected, 0.002);

    // alpha should be anti-symmetric
    double alpha_inv = ge.CalcAlpha(2, 1, 5.0, 1.0);
    EXPECT_NEAR(alpha_inv, -expected, 0.002);  //
  }
}

TEST(GroupedEvo, RungeKutta) {
  using vd_t = std::vector<double>;
  // dot{x} = -x
  {
    vd_t init = {1.0};
    std::function<void(const vd_t &, vd_t &)> calc_x_dot = [](const std::vector<double> &x,
                                                              std::vector<double> &x_dot) {
      x_dot[0] = -x[0];
    };
    auto x = GroupedEvo::SolveByRungeKutta(calc_x_dot, init, 0.01, 100, 0);
    // x(t) = 1.0 * exp(-t)
    EXPECT_NEAR(x[0], std::exp(-1.0), 0.002);
    x = GroupedEvo::SolveByRungeKutta(calc_x_dot, x, 0.01, 100, 0);
    EXPECT_NEAR(x[0], std::exp(-2.0), 0.002);
  }

  // Lotka-Volterra competition
  // dot{x} = x(3-x-2y)
  // dot{y} = y(2-x-y)
  // stable fixed points: (0,2), (3,0)
  {
    vd_t init = {1.0, 1.5};
    std::function<void(const vd_t &, vd_t &)> calc_x_dot = [](const std::vector<double> &x,
                                                              std::vector<double> &x_dot) {
      x_dot[0] = x[0] * (3.0 - x[0] - 2.0 * x[1]);
      x_dot[1] = x[1] * (2.0 - x[0] - x[1]);
    };
    auto x = GroupedEvo::SolveByRungeKutta(calc_x_dot, init, 0.1, 300, 0);
    EXPECT_NEAR(x[0], 0.0, 0.002);
    EXPECT_NEAR(x[1], 2.0, 0.002);

    init = {1.5, 1.0};
    x = GroupedEvo::SolveByRungeKutta(calc_x_dot, init, 0.1, 300, 0);
    EXPECT_NEAR(x[0], 3.0, 0.002);
    EXPECT_NEAR(x[1], 0.0, 0.002);
  }
}

TEST(GroupedEvo, SovleSampleSystem) {
  std::vector<size_t> norm_ids = {0, 1};
  std::vector<std::vector<double>> p_fix = {{0.0, 0.02}, {0.2, 0.0}};
  std::vector<double> self_coop_levels = {0.0, 1.0};
  GroupedEvo ge(norm_ids, p_fix, self_coop_levels);
  auto out = ge.TimeEvolutionODE(5.0, 1.0, 0.01, 1000, 1.0, std::cerr);
  IC(out);
  EXPECT_NEAR(out[0], 0.124, 0.002);
  EXPECT_NEAR(out[1], 0.876, 0.002);
}