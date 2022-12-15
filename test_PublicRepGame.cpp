#include <iostream>
#include <cassert>
#include <icecream.hpp>
#include "PublicRepGame.hpp"


template <typename T>
bool IsAllClose(T a, T b, double epsilon = 1.0e-6) {
  for (size_t i = 0; i < a.size(); i++) {
    if (std::abs(a[i] - b[i]) > epsilon) {
      return false;
    }
  }
  return true;
}

bool IsClose(double a, double b, double epsilon = 1.0e-6) {
  return std::abs(a - b) < epsilon;
}

void test_L8() {
  double mu_e = 0.001, mu_a_donor = 0.001, mu_a_recip = 0.001;
  std::vector<Norm>
      norms = {Norm::L1(), Norm::L2(), Norm::L3(), Norm::L4(), Norm::L5(), Norm::L6(), Norm::L7(), Norm::L8()};
  for (auto n : norms) {
    PublicRepGame g(mu_e, mu_a_donor, mu_a_recip, n);
    std::cerr << g.r_norm.Inspect();
    IC(g.h_star, g.pc_res_res);
    assert(g.h_star > 0.97);
    assert(g.pc_res_res > 0.96);

    ActionRule alld = ActionRule::ALLD();
    double H_alld = g.MutantEqReputation(alld);
    auto pcs_alld = g.MutantCooperationProbs(alld);
    IC(H_alld, pcs_alld);
    ActionRule allc = ActionRule::ALLC();
    double H_allc = g.MutantEqReputation(allc);
    auto pcs_allc = g.MutantCooperationProbs(allc);
    IC(H_allc, pcs_allc);
    assert(H_alld < 0.01);
    assert(H_allc > 0.99);

    auto br_alld = g.StableBenefitRangeAgainstMutant(alld);
    auto br_allc = g.StableBenefitRangeAgainstMutant(allc);
    IC(br_alld, br_allc);
    assert(br_alld[0] < 1.05);
    assert(br_allc[1] > 100);

    auto br = g.ESSBenefitRange();
    IC(br);
    assert(br[0] < 1.05);
    assert(br[1] > 100);
  }
  std::cout << "test_L8 passed!" << std::endl;
}

void test_ImageScoring() {
  double mu_e = 0.001, mu_a_donor = 0.001, mu_a_recip = 0.001;
  Norm is = Norm::ImageScoring();
  PublicRepGame g(mu_e, mu_a_donor, mu_a_recip, is);
  std::cerr << "IS: " << g.r_norm.Inspect();
  IC(g.h_star, g.pc_res_res);
  bool h = IsClose(g.h_star, 0.40, 0.01 );
  bool pc = IsClose( g.pc_res_res, 0.40, 0.01 );
  assert( h && pc );

  ActionRule alld = ActionRule::ALLD();
  ActionRule allc = ActionRule::ALLC();
  auto br_alld = g.StableBenefitRangeAgainstMutant(alld);
  auto br_allc = g.StableBenefitRangeAgainstMutant(allc);
  IC(br_alld, br_allc);

  auto br = g.ESSBenefitRange();
  IC(br);

  std::cout << "test_ImageScoring passed!" << std::endl;
}


int main() {
  test_L8();
  test_ImageScoring();

  return 0;
}