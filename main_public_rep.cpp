#include <iostream>
#include <cmath>
#include <vector>
#include <set>
#include <map>
#include <random>
#include <bitset>
#include <regex>
#include "Norm.hpp"
#include "PublicRepGame.hpp"


constexpr Reputation B = Reputation::B, G = Reputation::G;
constexpr Action C = Action::C, D = Action::D;

std::tuple<bool,std::array<double,2>,double> CheckCESS(const Norm& norm, double mu_a_recip = 1.0e-3) {
  const double mu_e = 1.0e-3, mu_a_donor = 1.0e-3;
  PublicRepGame game(mu_e, mu_a_donor, mu_a_recip, norm);
  auto brange = game.ESSBenefitRange();
  bool isCESS = game.pc_res_res > 0.98 && brange[0] < 10 && brange[0] + 0.001 < brange[1];
  return {isCESS, brange, game.h_star};
}

// comprehensively enumerate deterministic norms without R_2, and find CESSs
void FindLeadingEightSecondarySixteen() {
  size_t num = 0, num_passed = 0;

  for (int j = 0; j < 256; j++) {
    AssessmentRule R1 = AssessmentRule::MakeDeterministicRule(j);
    AssessmentRule R2 = AssessmentRule::KeepRecipient();

    for (int i = 0; i < 16; i++) {
      ActionRule P = ActionRule::MakeDeterministicRule(i);
      Norm norm(R1, R2, P);
      if ( norm.ID() < norm.SwapGB().ID() ) {
        continue;
      }

      auto [isCESS, brange, h_star] = CheckCESS(norm, 0.0);
      if (isCESS) {
        if (h_star < 0.5) {
          norm = norm.SwapGB();
        }
        std::string name = norm.GetName();
        // if name starts with "L"
        if (name[0] == 'L') {
          assert( std::abs(brange[0] - 1) < 0.01 );
        }
        else if (name[0] == 'S') {
          assert( std::abs(brange[0] - 2) < 0.01 );
        }
        else {
          throw std::runtime_error("unexpected norm");
        }
        assert( brange[1] > 10 );
        std::cerr << norm.Inspect();
        std::cerr << "brange: "  << brange[0] << " " << brange[1] << std::endl;
        num_passed++;
      }
      num++;
    }
  }

  assert(num_passed == 24);
  assert(num == 2080);
  IC(num, num_passed);
}

int IdentifyType(const Norm& norm) {
  std::string pr1_name = norm.PR1_Name();
  if (pr1_name.empty()) {
    return -1;
  }
  auto R1 = norm.Rd;
  auto R2 = norm.Rr;
  auto P = norm.P;

  if (R2.GProb(G,G,C) != 1) {
    return -1;
  }

  // if pr1_name matches regexp /^L[0-8]$/
  std::regex re("^L[1-8]$");
  if ( std::regex_match(pr1_name, re) ) {
    if (R2.GProb(G,B,D) == 0 && R2.GProb(B,G,C) == 1 ) {
      return 0;
    }
    else if (R2.GProb(G,B,D) == 1 && R2.GProb(B,G,C) == 0 ) {
      return 1;
    }
    else if (R2.GProb(G,B,D) == 1 && R2.GProb(B,G,C) == 1 ) {
      return 2;
    }
  }
  else if ( std::regex_match(pr1_name, std::regex("^S[0-9]+$")) ) {
    if (R2.GProb(G,B,D) == 0 && R2.GProb(B,G,D) == 1 ) {
      return 3;
    }
    else if (R2.GProb(G,B,D) == 1 && R2.GProb(B,G,D) == 0 ) {
      return 4;
    }
    else if (R2.GProb(G,B,D) == 1 && R2.GProb(B,G,D) == 1 ) {
      return 5;
    }
  }
  else if ( pr1_name == "L_prime" ) {
    if (R2.GProb(G,B,D) == 1 && R2.GProb(B,G,C) == 1 ) {
      return 6;
    }
  }
  else if ( pr1_name == "S_prime" ) {
    if (R2.GProb(G,B,D) == 1 && R2.GProb(B,G,D) == 1 ) {
      return 7;
    }
  }
  else if ( pr1_name == "S_prime_prime" ) {
    if (R2.GProb(G,B,D) == 1 && R2.GProb(B,G,D) == 1 ) {
      return 8;
    }
  }
  return -1;
}

void EnumerateAllDeterministicNorms() {
  size_t num = 0, num_passed = 0;
  std::map< int, std::vector<int>> cess_norms;

  for (int j = 0; j < 256; j++) {
    AssessmentRule R1 = AssessmentRule::MakeDeterministicRule(j);
    for (int k = 0; k < 256; k++) {
      AssessmentRule R2 = AssessmentRule::MakeDeterministicRule(k);
      for (int i = 0; i < 16; i++) {
        ActionRule P = ActionRule::MakeDeterministicRule(i);

        Norm norm(R1, R2, P);
        if ( norm.ID() < norm.SwapGB().ID() ) {
          continue;
        }

        auto [isCESS, brange, h_star] = CheckCESS(norm);
        if (isCESS) {
          if (h_star < 0.5) {
            norm = norm.SwapGB();
          }
          int type = IdentifyType(norm);
          // std::cerr << norm.Inspect(); IC(brange[0], brange[1], type);
          assert( type >= 0 );
          cess_norms[type].push_back(norm.ID());
          if      (type == 0) { assert(std::abs(brange[0] - 1) < 0.03); }
          else if (type == 1) { assert(std::abs(brange[0] - 2) < 0.03); }
          else if (type == 2) { assert(std::abs(brange[0] - 2) < 0.03); }
          else if (type == 3) { assert(std::abs(brange[0] - 2) < 0.03); }
          else if (type == 4) { assert(std::abs(brange[0] - 3) < 0.03); }
          else if (type == 5) { assert(std::abs(brange[0] - 3) < 0.03); }
          else if (type == 6) { assert(std::abs(brange[0] - 2) < 0.03); }
          else if (type == 7) { assert(std::abs(brange[0] - 2) < 0.03); }
          else if (type == 8) { assert(std::abs(brange[0] - 3) < 0.03); }
          num_passed++;
        }
        num++;
      }
    }
  }

  IC(num, num_passed);
  assert( num_passed == 2944 );

  for (int i = 0; i < 9; i++) {
    size_t s = cess_norms[i].size();
    // std::cerr << i << " " << s << std::endl;
    if (i >= 0 && i <= 2) {
      assert(s == 256);
    }
    else if (i >= 3 && i <= 5) {
      assert(s == 512);
    }
    else if (i == 6) {
      assert(s == 128);
    }
    else if (i >= 7) {
      assert(s == 256);
    }
  }
}

bool CloseEnough(double a, double b, double tol = 1.0e-2) {
  return std::abs(a - b) < tol;
}

std::array<double,2> AnalyticBenefitRange(const Norm& norm) {

  if (!norm.P.IsDeterministic()) {
    throw std::runtime_error("action rule is not deterministic");
  }

  auto R1 = [&norm](Reputation X, Reputation Y, Action A) -> double {
    return norm.Rd.GProb(X, Y, A);
  };
  auto R2 = [&norm](Reputation X, Reputation Y, Action A) -> double {
    return norm.Rr.GProb(X, Y, A);
  };
  auto P = [&norm](Reputation X, Reputation Y) -> double {
    return norm.P.CProb(X, Y);
  };

  // R_1(G,G,C) = 1 && R_2(G,G,C) = 1 &&
  // R_1(G,B) + R_2(G,B) + R_1(B,G) + R_2(B,G) > 2  => h* = 1
  Action bg = P(B, G) == 1.0 ? C : D;
  double r = R1(G,B,D) + R2(G,B,D)
      + R1(B,G,bg) + R2(B,G,bg);
  if (R1(G,G,C) == 1.0 &&
      R2(G,G,C) == 1.0 &&
      r > 2.0) {
    // h_star == 1 should be true
  }
  else {
    std::cerr << "R1(G,G,C) = " << R1(G,G,C) << std::endl;
    std::cerr << "R2(G,G,C) = " << R2(G,G,C) << std::endl;
    std::cerr << "r = " << r << std::endl;
    throw std::runtime_error("h_star is not 1");
  }

  // P(G,G) == 1
  if (P(G,G) != 1.0) {
    throw std::runtime_error("P(G,G) is not 1");
  }

  // P(G,B) == 0
  if (P(G,B) != 0.0) {
    throw std::runtime_error("P(G,B) == 0 is necessary");
  }

  std::array<double,2> b_range = {0.0, std::numeric_limits<double>::max()};

  // if P(B,G) = 1
  if (P(B,G) == 1.0) {

    // R_1(B,G,C) > R_1(B,G,D) &&
    // b/c > (R_1(B,G,C) + R_2(G,B) ) / (R_1(B,G,C) - R_1(B,G,D))
    if (R1(B,G,C) > R1(B,G,D) ) {
      double b_lower = (R1(B, G, C) + R2(G, B, D)) / (R1(B, G, C) - R1(B, G, D));
      std::cerr << "BG: b_lower = " << b_lower << std::endl;
      if (b_lower > b_range[0]) {
          b_range[0] = b_lower;
      }
    }
    else {
      b_range[0] = std::numeric_limits<double>::max();
      b_range[1] = 0.0;
      std::cerr << "BG: R_1(B,G,C) > R_1(B,G,D) is necessary" << std::endl;
      return b_range;
    }

    // R_1(G,G,D) < 1   &&
    // b/c > ( R_1(B,G,C) + R_2(G,B) ) / ( R_1(G,G,C) - R_1(G,G,D) )
    if (R1(G,G,D) < 1.0) {
      double b_lower = (R1(B, G, C) + R2(G, B, D)) / (R1(G, G, C) - R1(G, G, D));
      std::cerr << "GG: b_lower = " << b_lower << std::endl;
      if (b_lower > b_range[0]) {
          b_range[0] = b_lower;
      }
    }
    else {
      b_range[0] = std::numeric_limits<double>::max();
      b_range[1] = 0.0;
      std::cerr << "GG: R_1(G,G,D) < 1 is necessary" << std::endl;
      return b_range;
    }

    // R_1(G,B,C) <= R_1(G,B,D) ||
    // b/c < (R_1(B,G,C) + R_2(G,B) ) / ( R_1(G,B,C) - R_1(G,B,D) )
    if (R1(G,B,C) <= R1(G,B,D)) {
      double b_upper = std::numeric_limits<double>::infinity();
      std::cerr << "GB: b_upper = " << b_upper << std::endl;
    }
    else {
      double b_upper = (R1(B, G, C) + R2(G, B, D)) / (R1(G, B, C) - R1(G, B, D));
      std::cerr << "GB: b_upper = " << b_upper << std::endl;
      if (b_upper < b_range[1]) {
        b_range[1] = b_upper;
      }
    }

    // check P(B,B) is the optimal or not
    if ( P(B,B) == 1 ) {
      // R_1(B,B,C) > R_1(B,B,D)  &&
      // b/c > {R_1(B,G,C) + R_2(G,B)} / {R_1(B,B,C) - R_1(B,B,D)}
      if (R1(B,B,C) <= R1(B,B,D)) {
        b_range[0] = std::numeric_limits<double>::max();
        b_range[1] = 0.0;
        std::cerr << "BB: R_1(B,B,C) <= R_1(B,B,D) is necessary for P(B,B)=1 to be optimal" << std::endl;
        return b_range;
      }
      else {
        double b_lower = (R1(B, G, C) + R2(G, B, D)) / (R1(B, B, C) - R1(B, B, D));
        std::cerr << "BB: b_lower = " << b_lower << std::endl;
        if (b_lower > b_range[0]) {
          b_range[0] = b_lower;
        }
      }
    }
    else if ( P(B,B) == 0 ) {
      if (R1(B,B,C) <= R1(B,B,D)) { // always ESS
        std::cerr << "BB: always ESS" << std::endl;
      }
      else {
        double b_upper = (R1(B, G, C) + R2(G, B, D)) / (R1(B, B, C) - R1(B, B, D));
        std::cerr << "BB: b_upper = " << b_upper << std::endl;
        if (b_upper < b_range[1]) {
          b_range[1] = b_upper;
        }
      }
    }
  }
  else if (P(B,G) == 0) {

    // R_1(B,G,C) \leq R_1(B,G,D)  ||
    // b/c < {R_1(B,G,C)+R_2(G,B)} / {R_1(B,G,C)-R_1(B,G,D)}
    if (R1(B,G,C) <= R1(B,G,D) ) {
      std::cerr << "BG: always ESS" << std::endl;
    }
    else {
      double b_upper = (R1(B, G, C) + R2(G, B, D)) / (R1(B, G, C) - R1(B, G, D));
      std::cerr << "BG: b_upper = " << b_upper << std::endl;
      if (b_upper < b_range[1]) {
        b_range[1] = b_upper;
      }
    }

    // R_1(G,G,D) < 1 &&
    // b/c > 1 + {R_1(B,G,D)+R_2(G,B)} / {1-R_1(G,G,D)}
    if (R1(G,G,D) < 1.0) {
      double b_lower = 1.0 + (R1(B, G, D) + R2(G, B, D)) / (1.0 - R1(G, G, D));
      std::cerr << "GG: b_lower = " << b_lower << std::endl;
      if (b_lower > b_range[0]) {
        b_range[0] = b_lower;
      }
    }
    else {
      b_range[0] = std::numeric_limits<double>::max();
      b_range[1] = 0.0;
      std::cerr << "GG: R_1(G,G,D) < 1 is necessary" << std::endl;
    }

    // R_1(G,B,D) - R_1(G,B,C) \geq 0 ||
    // b/c < 1 + {R_1(B,G,D)+R_2(G,B)} / {R_1(G,B,C) - R_1(G,B,D)}
    if (R1(G,B,D) - R1(G,B,C) >= 0.0) {
      std::cerr << "GB: always ESS" << std::endl;
    }
    else {
      double b_upper = 1.0 + (R1(B, G, D) + R2(G, B, D)) / (R1(G, B, C) - R1(G, B, D));
      std::cerr << "GB: b_upper = " << b_upper << std::endl;
      if (b_upper < b_range[1]) {
        b_range[1] = b_upper;
      }
    }

    // check P(B,B) is the optimal or not
    if ( P(B,B) == 1 ) {
      // R_1(B,B,C) > R_1(B,B,D)  AND
      // b/c > 1 + {R_1(B,G,D) + R_2(G,B) } / {R_1(B,B,C) - R_1(B,B,D)}

      if (R1(B,B,C) <= R1(B,B,D)) {
        b_range[0] = std::numeric_limits<double>::max();
        b_range[1] = 0.0;
        std::cerr << "BB: R_1(B,B,C) <= R_1(B,B,D) is necessary for P(B,B)=1 to be optimal" << std::endl;
      }
      else {
        double b_lower = 1.0 + (R1(B, G, D) + R2(G, B, D)) / (R1(B, B, C) - R1(B, B, D));
        std::cerr << "BB: b_lower = " << b_lower << std::endl;
        if (b_lower > b_range[0]) {
          b_range[0] = b_lower;
        }
      }
    }
    else if ( P(B,B) == 0 ) {
      if (R1(B,B,C) <= R1(B,B,D)) { // always ESS
        std::cerr << "BB: always ESS" << std::endl;
      }
      else {
        double b_upper = 1.0 + (R1(B, G, D) + R2(G, B, D)) / (R1(B, B, C) - R1(B, B, D));
        std::cerr << "BB: b_upper = " << b_upper << std::endl;
        if (b_upper < b_range[1]) {
          b_range[1] = b_upper;
        }
      }
    }
  } else {
    throw std::runtime_error("P(B,G) == 0 or 1 is necessary");
  }

  return b_range;
}

Norm MakeNormFromTable(std::array<Action,4> actions, std::array<double,8> r1, std::array<double,8> r2 = {1,1,0,0,1,1,0,0}) {
  ActionRule ar({
                    {{G,G}, actions[0]==C?1:0},
                    {{G,B}, actions[1]==C?1:0},
                    {{B,G}, actions[2]==C?1:0},
                    {{B,B}, actions[3]==C?1:0},
  });

  AssessmentRule Rd({
                        {{G,G,C}, r1[0]}, {{G,G,D}, r1[1]},
                        {{G,B,C}, r1[2]}, {{G,B,D}, r1[3]},
                        {{B,G,C}, r1[4]}, {{B,G,D}, r1[5]},
                        {{B,B,C}, r1[6]}, {{B,B,D}, r1[7]},
  });
  AssessmentRule Rr({
                        {{G,G,C}, r2[0]}, {{G,G,D}, r2[1]},
                        {{G,B,C}, r2[2]}, {{G,B,D}, r2[3]},
                        {{B,G,C}, r2[4]}, {{B,G,D}, r2[5]},
                        {{B,B,C}, r2[6]}, {{B,B,D}, r2[7]},
                    });
  return Norm(Rd, Rr, ar);
}

void StochasticVariantLeadingEight() {
  {
    double p3 = 0.3;
    double p1 = 1.0 - p3;
    double p2 = 1.5 - p3;
    Norm norm = MakeNormFromTable({C, D, C, C}, {1, p1, 0, p2, p3, 0, 1, 0}, {1, 1, 0, 0, 1, 1, 0, 0});

    auto brange1 = AnalyticBenefitRange(norm);
    auto brange2 = PublicRepGame(1.0e-6, 1.0e-6, 1.0e-6, norm).ESSBenefitRange();
    IC(brange1, brange2);

    assert( std::abs(brange1[0] - 1) < 0.03 );
    assert( std::abs(brange2[0] - 1) < 0.03 );
    assert( brange1[1] > 10 && brange2[1] > 10 );
  }

  {
    double p3 = 0.2;
    double p1 = 0.0;
    double p2 = 1.0;
    Norm norm = MakeNormFromTable({C, D, D, D}, {1, p1, 0, p2, 0, p3, 0, 0}, {1, 1, 0, 0, 1, 1, 0, 0});
    auto brange1 = AnalyticBenefitRange(norm);
    auto brange2 = PublicRepGame(1.0e-6, 1.0e-6, 1.0e-6, norm).ESSBenefitRange();
    IC(brange1, brange2);

    assert( std::abs(brange1[0] - 1.2) < 0.03 );
    assert( std::abs(brange2[0] - 1.2) < 0.03 );
    assert( brange1[1] > 10 && brange2[1] > 10 );
  }
}

bool CompareAnalyticNumericalBranges(const Norm& norm) {
  std::cerr << norm.Inspect();
  auto brange1 = AnalyticBenefitRange(norm);
  auto brange2 = PublicRepGame(1.0e-6, 1.0e-6, 1.0e-6, norm).ESSBenefitRange();
  IC(brange1, brange2);

  if (brange1[0] >= brange1[1]) {  // no ESS
    return brange2[0] >= brange2[1];  // assert brange2 has no ESS too
  }

  // we don't care about unrealistically large b/c
  // it is susceptible to the numerical error
  const double large_th = 10.0;
  if (brange1[0] > large_th && brange2[0] > large_th) {
    // both have no ESS
    return true;
  }
  else if (brange1[0] > large_th && brange2[0] < large_th) {
    return false;
  }
  else if (brange1[0] < large_th && brange2[0] > large_th) {
    return false;
  }

  if (brange1[1] > large_th) { brange1[1] = large_th; }
  if (brange2[1] > large_th) { brange2[1] = large_th; }

  double th = 3.0e-2;
  return ( std::abs(brange1[0] - brange2[0]) < brange1[0]*th && std::abs(brange1[1] - brange2[1]) < brange1[1]*th );
}

void RandomCheckStochasticNorms() {
  size_t num_norms = 100'000;
  std::mt19937_64 rng(123456789);
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  auto r01 = [&rng,&dist]() -> double {
    return dist(rng);
  };

  std::vector<Norm> unpassed;
  for (size_t i = 0; i < num_norms; i++) {
    std::cerr << std::endl << "i: " << i << std::endl;
    Action P_bg = (rng() % 2 == 0) ? C : D;
    Action P_bb = (rng() % 2 == 0) ? C : D;
    std::array<Action,4> actions = {C,D,P_bg,P_bb};

    std::array<double,8> r1 = {
        1.0,   r01(), // GG
        r01(), r01(), // GB
        r01(), r01(), // BG
        r01(), r01()  // BB
    };
    std::array<double,8> r2 = {
        1.0,   r01(), // GG
        r01(), r01(), // GB
        r01(), r01(), // BG
        r01(), r01()  // BB
    };

    // R_1(G,B,D) + R_2(G,B,D) + R_1(B,G,C) + R_2(B,G,C) > 2 is necessary
    double recov = (P_bg==C) ? (r1[3] + r2[3] + r1[4] + r2[4]) : (r1[3] + r2[3] + r1[5] + r2[5]);
    if (recov <= 2.0) { std::cerr << "does not satisfy recov > 2. continue." << std::endl; continue; }

    Norm norm = MakeNormFromTable(actions, r1, r2);
    bool b = CompareAnalyticNumericalBranges(norm);
    if (!b) {
      unpassed.push_back(norm);
      // throw std::runtime_error("Failed");
    }

  }

  for (auto& norm : unpassed) {
    bool b = CompareAnalyticNumericalBranges(norm);
    std::cerr << "-----------------------------------------------------------\n";
    std::cerr << PublicRepGame(1.0e-6, 1.0e-6, 1.0e-6, norm).Inspect();
    std::cerr << "-----------------------------------------------------------\n\n";
  }
  std::cerr << "Number of unpassed: " << unpassed.size() << std::endl;
  assert(unpassed.size() == 0);

}

void RandomCheckSecondOrderStochasticNorms() {
  size_t num_norms = 10'000;
  std::mt19937_64 rng(9876543210);
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  auto r01 = [&rng,&dist]() -> double {
    return dist(rng);
  };

  for (size_t i = 0; i < num_norms; i++) {
    std::array<Action, 4> actions = {C, D, C, D};

    double r1gd = r01();  // R_1(-,G,D) < 1
    double r1bc = r01();
    double r1bd = r01();
    double r2gd = r01();
    double r2bc = r01();
    double r2bd = r01();
    std::array<double, 8> r1 = {
        1.0, r1gd, // GG
        r1bc, r1bd, // GB
        1.0, r1gd, // BG
        r1bc, r1bd  // BB
    };
    std::array<double, 8> r2 = {
        1.0, r2gd, // GG
        r2bc, r2bd, // GB
        1.0, r2gd, // BG
        r2bc, r2bd  // BB
    };

    // R_1(-,B,D) + R_2(-,B,D) > 0 is necessary
    double recov = r1[3] + r2[3];
    if (recov <= 1.0e-6) {
      std::cerr << "does not satisfy recov > 0. continue." << std::endl;
      continue;
    }
    if (1.0 - r1gd <= 1.0e-6) {
      std::cerr << "does not satisfy R1(-,G,D) < 1. continue." << std::endl;
      continue;
    }

    Norm norm = MakeNormFromTable(actions, r1, r2);
    double bc_min = (1.0 + r2bd) / (1.0 - r1gd);
    double bc_max = std::numeric_limits<double>::infinity();
    if (r1bc > r1bd) {
      bc_max = (1.0 + r2bd) / (r1bc - r1bd);
    }
    auto brange = AnalyticBenefitRange(norm);
    bool b = CompareAnalyticNumericalBranges(norm);

    if ( bc_min < 10 ) {
      assert(brange[0] < 10);
    }
    else {
      assert( std::abs(brange[0] - bc_min) < 1.0e-6 );
    }
    if ( bc_max > 10 ) {
      assert(brange[1] > 10);
    }
    else {
      assert( std::abs(brange[1] - bc_max) < 1.0e-6 );
    }

    IC(bc_min, bc_max, brange);
  }
}
std::vector<Norm> CESS_deterministic_norms(int c) {
  constexpr Reputation G = Reputation::G, B = Reputation::B;
  constexpr Action C = Action::C, D = Action::D;
  std::vector<Norm> ans;
  auto leading_eight = std::vector{Norm::L1(), Norm::L2(), Norm::L3(), Norm::L4(), Norm::L5(), Norm::L6(), Norm::L7(), Norm::L8()};
  auto secondary_sixteen = std::vector<Norm>();
  for (int i = 1; i <= 16; ++i) {
    secondary_sixteen.push_back(Norm::SecondarySixteen(i));
  }

  if (c == 0) {
    // leading eight + R_2(G,G,C) = 1, R_2(G,B,D) = 0, R_2(B,G,C) = 1
    for (auto &l : leading_eight) {
      for (int i = 0; i < 256; ++i) {
        AssessmentRule Rr = AssessmentRule::MakeDeterministicRule(i);
        if (Rr.GProb(G, G, C) == 1 && Rr.GProb(G, B, D) == 0 && Rr.GProb(B, G, C) == 1) {
          ans.push_back(Norm(l.Rd, Rr, l.P));
        }
      }
    }
    assert(ans.size() == 256);
  }
  else if (c == 1) {
    // leading eight + R_2(G,G,C) = 1, R_2(G,B,D) = 1, R_2(B,G,C) = 0
    for (auto &l : leading_eight) {
      for (int i = 0; i < 256; ++i) {
        AssessmentRule Rr = AssessmentRule::MakeDeterministicRule(i);
        if (Rr.GProb(G, G, C) == 1 && Rr.GProb(G, B, D) == 1 && Rr.GProb(B, G, C) == 0) {
          ans.push_back(Norm(l.Rd, Rr, l.P));
        }
      }
    }
  }
  else if (c == 2) {
    // leading eight + R_2(G,G,C) = 1, R_2(G,B,D) = 1, R_2(B,G,C) = 1
    for (auto &l : leading_eight) {
      for (int i = 0; i < 256; ++i) {
        AssessmentRule Rr = AssessmentRule::MakeDeterministicRule(i);
        if (Rr.GProb(G, G, C) == 1 && Rr.GProb(G, B, D) == 1 && Rr.GProb(B, G, C) == 1) {
          ans.push_back(Norm(l.Rd, Rr, l.P));
        }
      }
    }
  }
  else if (c == 3) {
    // secondary sixteen + R_2(G,G,C) = 1, R_2(G,B,D) = 0, R_2(B,G,D) = 1
    for (auto &s : secondary_sixteen) {
      for (int i = 0; i < 256; ++i) {
        AssessmentRule Rr = AssessmentRule::MakeDeterministicRule(i);
        if (Rr.GProb(G, G, C) == 1 && Rr.GProb(G, B, D) == 0 && Rr.GProb(B, G, D) == 1) {
          ans.push_back(Norm(s.Rd, Rr, s.P));
        }
      }
    }
  }
  else if (c == 4) {
    // secondary sixteen + R_2(G,G,C) = 1, R_2(G,B,D) = 1, R_2(B,G,D) = 0
    for (auto &s : secondary_sixteen) {
      for (int i = 0; i < 256; ++i) {
        AssessmentRule Rr = AssessmentRule::MakeDeterministicRule(i);
        if (Rr.GProb(G, G, C) == 1 && Rr.GProb(G, B, D) == 1 && Rr.GProb(B, G, D) == 0) {
          ans.push_back(Norm(s.Rd, Rr, s.P));
        }
      }
    }
  }
  else if (c == 5) {
    // secondary sixteen + R_2(G,G,C) = 1, R_2(G,B,D) = 1, R_2(B,G,D) = 1
    for (auto &s : secondary_sixteen) {
      for (int i = 0; i < 256; ++i) {
        AssessmentRule Rr = AssessmentRule::MakeDeterministicRule(i);
        if (Rr.GProb(G, G, C) == 1 && Rr.GProb(G, B, D) == 1 && Rr.GProb(B, G, D) == 1) {
          ans.push_back(Norm(s.Rd, Rr, s.P));
        }
      }
    }
  }
  else if (c == 6) {
    // L' [leading eight with R_1(G,B,C)==0 + R_1(G,B,D) = 0] + R_2(G,G,C) = 1, R_2(G,B,D) = 1, R_2(B,G,C) = 1
    std::vector<Norm> Lprime;
    for (auto l : leading_eight) {
      if (l.Rd.GProb(G, B, C) == 0.0) {
        l.Rd.SetGProb(G, B, D, 0.0);
        Lprime.push_back(l);
      }
    }
    assert(Lprime.size() == 4);
    for (auto l_prime: Lprime) {
      for (int i = 0; i < 256; ++i) {
        AssessmentRule Rr = AssessmentRule::MakeDeterministicRule(i);
        if (Rr.GProb(G, G, C) == 1 && Rr.GProb(G, B, D) == 1 && Rr.GProb(B, G, C) == 1) {
          ans.push_back(Norm(l_prime.Rd, Rr, l_prime.P));
        }
      }
    }
  }
  else if (c == 7) {
    // S' [secondary sixteen with R_1(B,G,C)==0 + R_1(B,G,D) = 0] + R_2(G,G,C) = 1, R_2(G,B,D) = 1, R_2(B,G,D) = 1
    std::vector<Norm> Sprime;
    for (auto s : secondary_sixteen) {
      if (s.Rd.GProb(B, G, C) == 0.0) {
        s.Rd.SetGProb(B, G, D, 0.0);
        Sprime.push_back(s);
      }
    }
    assert(Sprime.size() == 8);
    for (auto s_prime: Sprime) {
      for (int i = 0; i < 256; ++i) {
        AssessmentRule Rr = AssessmentRule::MakeDeterministicRule(i);
        if (Rr.GProb(G, G, C) == 1 && Rr.GProb(G,B,D) == 1 && Rr.GProb(B, G, D) == 1) {
          ans.push_back(Norm(s_prime.Rd, Rr, s_prime.P));
        }
      }
    }
  }
  else if (c == 8) {
    // S'' [secondary sixteen with R_1(G,B,C)==0 + R_1(G,B,D) = 0] + R_2(G,G,C) = 1, R_2(G,B,D) = 1, R_2(B,G,D) = 1
    std::vector<Norm> Sprime;
    for (auto s : secondary_sixteen) {
      if (s.Rd.GProb(G,B, C) == 0.0) {
        s.Rd.SetGProb(G,B, D, 0.0);
        Sprime.push_back(s);
      }
    }
    assert(Sprime.size() == 8);
    for (auto s_prime : Sprime) {
      for (int i = 0; i < 256; ++i) {
        AssessmentRule Rr = AssessmentRule::MakeDeterministicRule(i);
        if (Rr.GProb(G, G, C) == 1 && Rr.GProb(G, B, D) == 1 && Rr.GProb(B, G, D) == 1) {
          ans.push_back(Norm(s_prime.Rd, Rr, s_prime.P));
        }
      }
    }
  }
  else {
    throw std::runtime_error("invalid c");
  }
  return ans;
}

void check_CESS_deterministic_norms() {
  auto close = [](double a, double b) {
    return std::abs(a - b) < 0.05;
  };

  auto assert_CESS = [&close](int cess_class, int num_norms, double b_lower, int num_second_oder) {
    auto norms = CESS_deterministic_norms(cess_class);
    assert(norms.size() == num_norms);
    std::vector<Norm> second_order_norms;
    for (auto& norm : norms) {
      auto [isCESS, brange, h_star] = CheckCESS(norm);
      if (!isCESS) {
        std::cerr << "Failed CESS class " << cess_class << " with " << num_norms << " norms" << std::endl;
        IC(isCESS, brange, h_star, norm.Inspect());
      }
      assert(isCESS);
      assert(close(brange[0], b_lower));
      assert(brange[1] > 1.0e2);
      assert(close(h_star, 1.0));
      int type = IdentifyType(norm);
      assert(type == cess_class);
      if ( norm.IsSecondOrder() ) {
        second_order_norms.push_back(norm);
      }
    }
    assert(second_order_norms.size() == num_second_oder);
    std::cerr << "Passed CESS class " << cess_class << " with " << num_norms << " norms" << std::endl;
  };

  assert_CESS(0, 256, 1.0, 8);
  assert_CESS(1, 256, 2.0, 0);
  assert_CESS(2, 256, 2.0, 8);
  assert_CESS(3, 512, 2.0, 0);
  assert_CESS(4, 512, 3.0, 0);
  assert_CESS(5, 512, 3.0, 0);
  assert_CESS(6, 128, 2.0, 4);
  assert_CESS(7, 256, 2.0, 0);
  assert_CESS(8, 256, 3.0, 0);

}

void CESS_coop_classification() {
  // conduct classification of the norms based on the cooperation level at mu=1e-3

  using key_type = std::pair<int,double>; // (b_c_lower, error_sensitivity)
  std::map<key_type, std::vector<Norm> > norm_map;

  auto check_norms = [&norm_map](int norm_class, int exp_b_lower, double exp_error_sens1, double exp_error_sens2) {
    double mu = 1.0e-4;
    auto norms = CESS_deterministic_norms(norm_class);
    for (auto &norm : norms) {
      PublicRepGame game(mu, mu, mu, norm);
      double bc_lower = game.ESSBenefitRange()[0];
      double error_sensitivity = (1.0 - game.pc_res_res) / mu;
      error_sensitivity = std::round(error_sensitivity * 10.0) / 10.0;
      auto key = std::make_pair((int)std::round(bc_lower), error_sensitivity);
      norm_map[key].push_back(norm);

      assert( key.first == exp_b_lower );
      if (norm.Rr.GProb(G, G, D) == 1.0) {
        assert( key.second == exp_error_sens1 );
      }
      else if (norm.Rr.GProb(G, G, D) == 0.0) {
        assert( key.second == exp_error_sens2 );
      }
      else {
        throw std::runtime_error("unexpected cooperation level");
      }
    }
  };

  check_norms(0, 1, 4.0, 5.0);
  check_norms(1, 2, 4.0, 5.0);
  check_norms(2, 2, 2.5, 3.0);
  check_norms(3, 2, 7.0, 9.0);
  check_norms(4, 3, 7.0, 9.0);
  check_norms(5, 3, 4.0, 5.0);
  check_norms(6, 2, 4.0, 5.0);
  check_norms(7, 2, 7.0, 9.0);
  check_norms(8, 3, 7.0, 9.0);
  // map the values in norm_map to its size
  std::map<key_type, int> size_map;
  for (auto& kv : norm_map) {
    size_map[kv.first] = kv.second.size();
  }
  IC(size_map);
}

int main() {

  FindLeadingEightSecondarySixteen();
  std::cerr << "FindLeadingEightSecondarySixteen() done" << std::endl;
  EnumerateAllDeterministicNorms();
  std::cerr << "EnumerateAllDeterministicNorms() done" << std::endl;
  StochasticVariantLeadingEight();
  std::cerr << "StochasticVariantLeadingEight() done" << std::endl;
  RandomCheckStochasticNorms();
  std::cerr << "RandomCheckStochasticNorms() done" << std::endl;
  RandomCheckSecondOrderStochasticNorms();
  std::cerr << "RandomCheckSecondOrderStochasticNorms() done" << std::endl;
  check_CESS_deterministic_norms();
  std::cerr << "check_CESS_deterministic_norms() done" << std::endl;
  CESS_coop_classification();
  std::cerr << "CESS_coop_classification() done" << std::endl;

  return 0;
}
