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
