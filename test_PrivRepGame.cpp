#include <iostream>
#include <fstream>
#include <chrono>
#include <regex>
#include <icecream.hpp>
#include <nlohmann/json.hpp>
#include "PrivRepGame.hpp"

#include <gtest/gtest.h>

void test_SelfCooperationLevel(const Norm& norm, double expected_c_level, double expected_good_rep) {
  PrivateRepGame priv_game( {{norm, 50}}, 123456789ull);
  priv_game.Update(1e4, 0.9, 0.05, 0.0, false);
  priv_game.ResetCounts();
  priv_game.Update(1e4, 0.9, 0.05, 0.0, true);
  // IC( priv_game.NormCooperationLevels(), priv_game.NormAverageReputation() );
  EXPECT_NEAR( priv_game.SystemWideCooperationLevel(), expected_c_level, 0.02);
  EXPECT_NEAR( priv_game.NormCooperationLevels()[0][0], expected_c_level, 0.02);
  EXPECT_NEAR( priv_game.NormAverageReputation()[0][0], expected_good_rep, 0.02);
}

TEST(SelfCooperationLevel, RandomNorm) {
  test_SelfCooperationLevel(Norm::Random(), 0.5, 0.5);
}

TEST(SelfCooperationLevel, LeadingEight) {
  test_SelfCooperationLevel(Norm::L1(), 0.90, 0.90);
  test_SelfCooperationLevel(Norm::L2(), 0.66, 0.65);
  test_SelfCooperationLevel(Norm::L3(), 0.90, 0.90);
  test_SelfCooperationLevel(Norm::L4(), 0.90, 0.90);
  test_SelfCooperationLevel(Norm::L5(), 0.70, 0.70);
  test_SelfCooperationLevel(Norm::L6(), 0.50, 0.50);
  test_SelfCooperationLevel(Norm::L7(), 0.88, 0.88);
  test_SelfCooperationLevel(Norm::L8(), 0.0, 0.0);
}

TEST(PrivateRepGame, RandomNonIdenticalPermutations) {
  PrivateRepGame priv_game({{Norm::Random(), 3}}, 123456789ull);
  for (int t = 0; t < 100; t++) {
    std::mt19937_64 rng(123456789ull + t);
    auto [perm1, perm2] = priv_game.RandomNonIdenticalPermutations(2, rng);
    EXPECT_EQ(perm1.size(), 2);
    EXPECT_EQ(perm2.size(), 2);
    std::set<size_t> s1, s2;
    for (int i = 0; i < perm1.size(); i++) {
      EXPECT_NE(perm1[i], perm2[i]);
      s1.insert(perm1[i]);
      s2.insert(perm2[i]);
    }
    EXPECT_EQ(s1.size(), 2);
    EXPECT_EQ(s2.size(), 2);
  }
  for (int t = 0; t < 100; t++) {
    std::mt19937_64 rng(123456789ull + t);
    auto [perm1, perm2] = priv_game.RandomNonIdenticalPermutations(3, rng);
    EXPECT_EQ(perm1.size(), 3);
    EXPECT_EQ(perm2.size(), 3);
    std::set<size_t> s1, s2;
    for (int i = 0; i < perm1.size(); i++) {
      EXPECT_NE(perm1[i], perm2[i]);
      s1.insert(perm1[i]);
      s2.insert(perm2[i]);
    }
    EXPECT_EQ(s1.size(), 3);
    EXPECT_EQ(s2.size(), 3);
  }
  for (int t = 0; t < 100; t++) {
    std::mt19937_64 rng(123456789ull + t);
    auto [perm1, perm2] = priv_game.RandomNonIdenticalPermutations(50, rng);
    EXPECT_EQ(perm1.size(), 50);
    EXPECT_EQ(perm2.size(), 50);
    std::set<size_t> s1, s2;
    for (int i = 0; i < perm1.size(); i++) {
      EXPECT_NE(perm1[i], perm2[i]);
      s1.insert(perm1[i]);
      s2.insert(perm2[i]);
    }
    EXPECT_EQ(s1.size(), 50);
    EXPECT_EQ(s2.size(), 50);
  }
}

TEST(EvolPrivRepGame, L1_AllC_AllD) {
  Norm norm = Norm::L1();
  auto params = EvolPrivRepGame::SimulationParameters::Default();
  params.n_init = 1e5;
  params.n_steps = 1e5;

  EvolPrivRepGame evol(params);
  auto rhos = evol.FixationProbabilities({norm, Norm::AllC(), Norm::AllD()}, 5.0, 1.0);
  // rhos: [
  // [0, 0.0976933, 4.47544e-29],
  // [0.0114276, 0, 0.668377],
  // [0.0444453, 2.16425e-24, 0]
  // ]

  EXPECT_NEAR(rhos[0][1], 0.097, 0.02);
  EXPECT_NEAR(rhos[0][2], 0.000, 0.02);
  EXPECT_NEAR(rhos[1][0], 0.012, 0.02);
  EXPECT_NEAR(rhos[2][0], 0.043, 0.02);
  EXPECT_NEAR(rhos[1][2], 0.668, 0.02);
  EXPECT_NEAR(rhos[2][1], 0.000, 0.02);

  auto eq = evol.EquilibriumPopulationLowMut(rhos);
  // eq: [0.302589, 0.0434844, 0.653927]
  EXPECT_NEAR(eq[0], 0.30, 0.02);
  EXPECT_NEAR(eq[1], 0.04, 0.02);
  EXPECT_NEAR(eq[2], 0.66, 0.02);
}

// TEST(EvolPrivRepGame, FixationProbability) {
//   EvolPrivRepGame::SimulationParameters params;
//   params.n_init = 1e5;
//   params.n_steps = 1e5;
//
//   EvolPrivRepGame evol(50, {Norm::L1(), Norm::AllC(), Norm::AllD()}, params);
//
// }

TEST(EvolPrivRepGameAllCAllD, L1) {
  auto params = EvolPrivRepGame::SimulationParameters::Default();
  params.n_init = 1e5;
  params.n_steps = 1e5;

  EvolPrivRepGameAllCAllD evol(params, 5.0, 1.0);

  auto selfc_rho_eq = evol.EquilibriumCoopLevelAllCAllD(Norm::L1());
  double self_cooperation_level = std::get<0>(selfc_rho_eq);
  auto rhos = std::get<1>(selfc_rho_eq);
  auto eq = std::get<2>(selfc_rho_eq);

  // IC(self_cooperation_level, rhos, eq);
  // ic| self_cooperation_level: 0.90224
  // ic| rhos: [
  // [0, 0.0933583, 4.4365e-29],
  // [0.0120976, 0, 0.670243],
  // [0.0453502, 2.35795e-24, 0]
  // ]
  // ic| eq: [0.316563, 0.0433123, 0.640125]

  EXPECT_NEAR(self_cooperation_level, 0.90, 0.02);
  EXPECT_NEAR(rhos[0][1], 0.097, 0.02);
  EXPECT_NEAR(rhos[0][2], 0.000, 0.02);
  EXPECT_NEAR(rhos[1][0], 0.012, 0.02);
  EXPECT_NEAR(rhos[2][0], 0.043, 0.02);
  EXPECT_NEAR(eq[0], 0.30, 0.02);
  EXPECT_NEAR(eq[1], 0.04, 0.02);
  EXPECT_NEAR(eq[2], 0.66, 0.02);
}

TEST(EvolPrivRepGameFiniteMutationRateAllCAllD, L1) {
  auto params = EvolPrivRepGame::SimulationParameters::Default();
  params.n_init = 1e4;
  params.n_steps = 1e4;

  Norm L1 = Norm::L1();
  EvolPrivRepGameFiniteMutationRateAllCAllD evol(L1, params);
  double mu = 1.0e-2;
  double benefit = 5.0;
  auto result = evol.CalculateEquilibrium(benefit, 1.0, mu);
  EXPECT_NEAR(result.OverallCooperationLevel(), 0.74, 0.02);
  EXPECT_NEAR(result.OverallAbundances()[0], 0.47, 0.02);
  EXPECT_NEAR(result.OverallAbundances()[1], 0.32, 0.02);
  EXPECT_NEAR(result.OverallAbundances()[2], 0.21, 0.02);
}
