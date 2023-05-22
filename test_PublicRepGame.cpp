#include <iostream>
#include <regex>
#include <icecream.hpp>
#include "PublicRepGame.hpp"

#include <gtest/gtest.h>

TEST(PublicRepGame, L8) {
  double mu_e = 0.001, mu_a_donor = 0.001, mu_a_recip = 0.001;
  std::vector<Norm>
      norms = {Norm::L1(), Norm::L2(), Norm::L3(), Norm::L4(), Norm::L5(), Norm::L6(), Norm::L7(), Norm::L8()};
  for (auto n : norms) {
    PublicRepGame g(mu_e, mu_a_donor, mu_a_recip, n);
    // std::cerr << g.r_norm.Inspect();
    // IC(g.h_star, g.pc_res_res);
    EXPECT_GT(g.h_star, 0.97);
    EXPECT_GT(g.pc_res_res, 0.96);

    ActionRule alld = ActionRule::ALLD();
    double H_alld = g.MutantEqReputation(alld);
    EXPECT_LT(H_alld, 0.01);
    ActionRule allc = ActionRule::ALLC();
    double H_allc = g.MutantEqReputation(allc);
    EXPECT_GT(H_allc, 0.99);

    auto br_alld = g.StableBenefitRangeAgainstMutant(alld);
    auto br_allc = g.StableBenefitRangeAgainstMutant(allc);
    // IC(br_alld, br_allc);
    // ic| br_alld: [1.00201, 1.79769e+308]
    //     br_allc: [-1.00803, 1.79769e+308]
    EXPECT_LT(br_alld[0], 1.05);
    EXPECT_GT(br_allc[1], 100);

    auto br = g.ESSBenefitRange();
    // IC(br);
    // ic| br: [1.00201, 1.79769e+308]
    EXPECT_LT(br[0], 1.05);
    EXPECT_GT(br[1], 100);
  }
}

TEST(PublicRepGame, ImageScoring) {
  double mu_e = 0.001, mu_a_donor = 0.001, mu_a_recip = 0.001;
  Norm is = Norm::ImageScoring();
  PublicRepGame g(mu_e, mu_a_donor, mu_a_recip, is);
  // std::cerr << "IS: " << g.r_norm.Inspect();
  EXPECT_NEAR(g.h_star, 0.40, 0.01 );
  EXPECT_NEAR( g.pc_res_res, 0.40, 0.01 );

  ActionRule alld = ActionRule::ALLD();
  ActionRule allc = ActionRule::ALLC();
  auto br_alld = g.StableBenefitRangeAgainstMutant(alld);
  auto br_allc = g.StableBenefitRangeAgainstMutant(allc);
  // ic| br_alld: [1.00501, 1.79769e+308], br_allc: [0, 1.00501]
  // stable against allD but not against AllC
  EXPECT_NEAR(br_alld[0], 1.005, 0.01);
  EXPECT_GT(br_alld[1], 100);
  EXPECT_NEAR(br_allc[1], 1, 0.1);

  auto br = g.ESSBenefitRange();
  EXPECT_LT(br[1]-br[0], 0.1);
}

void PrintESSRange(const Norm& norm) {
  double mu_e = 0.001, mu_a_donor = 0.001, mu_a_recip = 0.001;
  PublicRepGame g(mu_e, mu_a_donor, mu_a_recip, norm);
  std::cerr << g.r_norm.Inspect();
  std::cerr << "h*: " << g.h_star << ", pc_res_res: " << g.pc_res_res << std::endl;

  ActionRule alld = ActionRule::ALLD();
  auto br_alld = g.StableBenefitRangeAgainstMutant(alld);
  std::cerr << "stable benefit range against ALLD: " << br_alld[0] << ", " << br_alld[1] << std::endl;
  ActionRule allc = ActionRule::ALLC();
  auto br_allc = g.StableBenefitRangeAgainstMutant(allc);
  std::cerr << "stable benefit range against ALLC: " << br_alld[0] << ", " << br_alld[1] << std::endl;

  auto br = g.ESSBenefitRange();
  std::cerr << "ESS b_range: " << br[0] << ", " << br[1] << std::endl;
}

/*
int main(int argc, char* argv[]) {

  std::vector<std::string> args;
  bool swap_gb = false;
  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == "-s") {
      swap_gb = true;
    }
    else {
      args.emplace_back(argv[i]);
    }
  }

  if (args.empty()) {
    test_L8();
    test_ImageScoring();
  }
  else if (args.size() == 1) {
    Norm n = Norm::ParseNormString(argv[1], swap_gb);
    PrintESSRange(n);
  }
  else {
    std::cerr << "Unsupported usage" << std::endl;
    return 1;
  }

  return 0;
}
 */