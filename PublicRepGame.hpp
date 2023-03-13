#ifndef GAME_HPP
#define GAME_HPP

#include <iostream>
#include <sstream>
#include <cmath>
#include <functional>
#include <iomanip>
#include <icecream.hpp>
#include "Norm.hpp"

class PublicRepGame {
public:
  PublicRepGame(double mu_e, double mu_a_donor, double mu_a_recip, const Norm& norm) :
  mu_e(mu_e), mu_a_donor(mu_a_donor), mu_a_recip(mu_a_recip), h_star(0.0), pc_res_res(0.0),
  norm(norm), r_norm(norm.RescaleWithError(mu_e, mu_a_donor, mu_a_recip)) {
    const_cast<double&>(h_star) = ResidentEqReputation();
    const_cast<double&>(pc_res_res) = h_star * h_star * r_norm.CProb(Reputation::G, Reputation::G)
        + h_star * (1.0-h_star) * r_norm.CProb(Reputation::G, Reputation::B)
        + (1.0-h_star) * h_star * r_norm.CProb(Reputation::B, Reputation::G)
        + (1.0-h_star) * (1.0-h_star) * r_norm.CProb(Reputation::B, Reputation::B);
  }
  std::string Inspect() const {
    std::stringstream ss;
    ss << "PublicRepGame:" << std::endl
       << "(mu_e, mu_a_donor, mu_a_recip): (" << mu_e << ", " << mu_a_donor << ", " << mu_a_recip << ")" << std::endl
       << "(Norm): " << norm.Inspect()
       << "(c_prob,h*): " << pc_res_res << ' ' << h_star << std::endl;
    for (int id = 0; id < 16; id++) {
      if (norm.P.ID() == id) continue;
      ActionRule mut = ActionRule::MakeDeterministicRule(id);
      auto pc = MutantCooperationProbs(mut);
      auto b_range = StableBenefitRangeAgainstMutant(mut);
      ss << "Against mutant: ";
      ss << ((mut.CProb(Reputation::G, Reputation::G) == 1.0) ? "C" : "D");
      ss << ((mut.CProb(Reputation::G, Reputation::B) == 1.0) ? "C" : "D");
      ss << ((mut.CProb(Reputation::B, Reputation::G) == 1.0) ? "C" : "D");
      ss << ((mut.CProb(Reputation::B, Reputation::B) == 1.0) ? "C" : "D");
      ss << " , benefit range: " << b_range[0] << ' ' << b_range[1] << " , (pc_res_mut,pc_mut_res): " << pc.first << ", " << pc.second << std::endl;
    }
    return ss.str();
  }
  const double mu_e, mu_a_donor, mu_a_recip;
  const double h_star; // equilibrium reputation of resident species
  const double pc_res_res;  // cooperation probability of resident species
  const Norm norm;
  const Norm r_norm;  // norm that takes into account error rates (mu_e, mu_a)
private:
  double P(Reputation X, Reputation Y) const { return r_norm.CProb(X, Y); }
  double R1(Reputation X, Reputation Y, Action A) const { return r_norm.GProbDonor(X, Y, A); }
  double R2(Reputation X, Reputation Y, Action A) const { return r_norm.GProbRecip(X, Y, A); }
public:
  double ResidentEqReputation() const {
    // dynamcis of h(t)
    // dot{h(t)} =
    //   h^2    { R1(G,G) + R2(G,G) }
    // + h(1-h) { R1(G,B) + R2(G,B) + R1(B,G) + R2(B,G) }
    // + (1-h)^2{ R1(B,B) + R2(B,B) }
    // -2h(t)

    constexpr Reputation B = Reputation::B, G = Reputation::G;
    constexpr Action C = Action::C, D = Action::D;

    double R1_BB = P(B, B) * R1(B, B, C) + (1.0 - P(B, B)) * R1(B, B, D);
    double R1_BG = P(B, G) * R1(B, G, C) + (1.0 - P(B, G)) * R1(B, G, D);
    double R1_GB = P(G, B) * R1(G, B, C) + (1.0 - P(G, B)) * R1(G, B, D);
    double R1_GG = P(G, G) * R1(G, G, C) + (1.0 - P(G, G)) * R1(G, G, D);

    double R2_BB = P(B, B) * R2(B, B, C) + (1.0 - P(B, B)) * R2(B, B, D);
    double R2_BG = P(B, G) * R2(B, G, C) + (1.0 - P(B, G)) * R2(B, G, D);
    double R2_GB = P(G, B) * R2(G, B, C) + (1.0 - P(G, B)) * R2(G, B, D);
    double R2_GG = P(G, G) * R2(G, G, C) + (1.0 - P(G, G)) * R2(G, G, D);

    // IC(R1_GG, R1_GB, R1_BG, R1_BB);
    // IC(R2_GG, R2_GB, R2_BG, R2_BB);

    // A = R_1(G,G)+R_2(G,G) -R_1(G,B)-R_2(G,B) -R_1(B,G)-R_2(B,G) + R_1(B,B)+R_2(B,B)
    double a = R1_GG + R2_GG - R1_GB - R2_GB - R1_BG - R2_BG + R1_BB + R2_BB;
    // B = R_1(G,B)+R_2(G,B) +R_1(B,G)+R_2(B,G) -2R_1(B,B) -2R_2(B,B) -2
    double b = R1_GB + R2_GB + R1_BG + R2_BG - 2.0 * R1_BB - 2.0 * R2_BB - 2.0;
    // C = R_1(B,B)+R_2(B,B)
    double c = R1_BB + R2_BB;

    // Although analytic expression is available, it is numerically unstable when A is small.
    // Use newton's method instead.
    const double tolerance = 2.0e-3;
    double h = (std::abs(a) > tolerance) ? ((-b - std::sqrt(b*b - 4.0*a*c)) / (2.0*a)) : -c / b;
    // solve f(h) = ah^2 + bh + c = 0 by Newton's method
    if (h > 1.0) h = 1.0;
    if (h < 0.0) h = 0.0;
    double f = a*h*h + b*h + c;
    double df = 2.0*a*h + b;
    while (std::abs(f) > 1.0e-12) {
      h -= f / df;
      f = a*h*h + b*h + c;
      df = 2.0*a*h + b;
    }
    // IC(h,f);
    return h;
  }

  std::array<double,2> ESSBenefitRange() const {  // range of b for which the norm is ESS against possible mutants
    double b_lower_bound = 1.0;
    double b_upper_bound = std::numeric_limits<double>::max();
    for (int id = 0; id < 16; id++) {
      if (norm.P.ID() == id) continue;
      ActionRule mut = ActionRule::MakeDeterministicRule(id);
      auto b_range = StableBenefitRangeAgainstMutant(mut);
      if (b_range[0] > b_lower_bound) b_lower_bound = b_range[0];
      if (b_range[1] < b_upper_bound) b_upper_bound = b_range[1];
    }
    return {b_lower_bound, b_upper_bound};
  }

  std::array<double,2> StableBenefitRangeAgainstMutant(const ActionRule& mut) const {
    double b_lower_bound = 1.0;
    double b_upper_bound = std::numeric_limits<double>::max();

    auto p = MutantCooperationProbs(mut);
    double pc_res_mut = p.first, pc_mut_res = p.second;

    if (pc_res_res > pc_res_mut) {
      // b/c > { p_c^{\rm res \to res} - p_c^{\rm mut \to res} / { p_c^{\rm res \to res} - p_c^{\rm res \to mut} }
      b_lower_bound = (pc_res_res - pc_mut_res) / (pc_res_res - pc_res_mut);
      b_upper_bound = std::numeric_limits<double>::max();
    }
    else if (pc_res_res < pc_res_mut) {
      // b/c < { p_c^{\rm res \to mut} - p_c^{\rm mut \to res} / { p_c^{\rm res \to res} - p_c^{\rm res \to mut} }
      b_upper_bound = (pc_res_res - pc_mut_res) / (pc_res_res - pc_res_mut);
      b_lower_bound = 0.0;
    }
    else {
      if (pc_res_res >= pc_mut_res) {
        // no ESS
        b_lower_bound = std::numeric_limits<double>::max();
        b_upper_bound = 0.0;
      }
      else {
        // ESS for all
        b_lower_bound = 1.0;
        b_upper_bound = std::numeric_limits<double>::max();
      }
    }
    return {b_lower_bound, b_upper_bound};
  }

  double MutantPayoff(const ActionRule& mut, double benefit) const {
    auto pcs = MutantCooperationProbs(mut);
    double pc_res_mut = pcs.first;
    double pc_mut_res = pcs.second;
    double pi_mut = pc_res_mut * benefit - pc_mut_res * 1.0;
    return pi_mut;
  }

  std::pair<double,double> MutantCooperationProbs(const ActionRule& mut) const {
    constexpr Reputation B = Reputation::B, G = Reputation::G;
    constexpr Action C = Action::C, D = Action::D;
    double H_star = MutantEqReputation(mut);
    const ActionRule Pmut = mut.RescaleWithError(mu_e);
    double pc_res_mut = h_star * H_star * P(G, G)
        + h_star * (1.0-H_star) * P(G, B)
        + (1.0-h_star) * H_star * P(B, G)
        + (1.0-h_star) * (1.0-H_star) * P(B, B);
    double pc_mut_res = H_star * h_star * Pmut.CProb(G, G)
        + H_star * (1.0-h_star) * Pmut.CProb(G, B)
        + (1.0-H_star) * h_star * Pmut.CProb(B, G)
        + (1.0-H_star) * (1.0-h_star) * Pmut.CProb(B, B);
    return std::make_pair(pc_res_mut, pc_mut_res);
  }

  double MutantEqReputation(const ActionRule& mut) const {
    // H^{\ast} =
    //   h^{\ast} [ R_1(B,G,P') + R_2(G,B,P) ] + (1-h^{\ast}) [ R_1(B,B,P') + R_2(B,B,P)
    //   / { 2 - h^{\ast} [R_1(G,G,P') + R_2(G,G,P) - R_1(B,G,P') - R_2(G,B,P) ] - (1-h^{\ast}) [ R_1(G,B,P') + R_2(B,G,P) - R_1(B,B,P') - R_2(B,B,P)] }.
    constexpr Reputation B = Reputation::B, G = Reputation::G;
    constexpr Action C = Action::C, D = Action::D;
    const ActionRule Pmut = mut.RescaleWithError(mu_e);
    auto Pm = [Pmut](Reputation X, Reputation Y)->double { return Pmut.CProb(X, Y); };

    double r1_BBPmut = Pm(B, B) * R1(B, B, C) + (1.0 - Pm(B, B)) * R1(B, B, D);
    double r1_BGPmut = Pm(B, G) * R1(B, G, C) + (1.0 - Pm(B, G)) * R1(B, G, D);
    double r1_GBPmut = Pm(G, B) * R1(G, B, C) + (1.0 - Pm(G, B)) * R1(G, B, D);
    double r1_GGPmut = Pm(G, G) * R1(G, G, C) + (1.0 - Pm(G, G)) * R1(G, G, D);

    double r2_BBPres = P(B, B) * R2(B, B, C) + (1.0 - P(B, B)) * R2(B, B, D);
    double r2_BGPres = P(B, G) * R2(B, G, C) + (1.0 - P(B, G)) * R2(B, G, D);
    double r2_GBPres = P(G, B) * R2(G, B, C) + (1.0 - P(G, B)) * R2(G, B, D);
    double r2_GGPres = P(G, G) * R2(G, G, C) + (1.0 - P(G, G)) * R2(G, G, D);

    double num = h_star * (r1_BGPmut + r2_GBPres) + (1.0 - h_star) * (r1_BBPmut + r2_BBPres);
    double den = 2.0 - h_star * (r1_GGPmut + r2_GGPres - r1_BGPmut - r2_GBPres) - (1.0 - h_star) * (r1_GBPmut + r2_BGPres - r1_BBPmut - r2_BBPres);
    double H_star = num / den;
    double H_dot = ( h_star * H_star ) * (r1_GGPmut + r2_GGPres)
        + ( h_star * (1.0 - H_star) ) * (r1_BGPmut + r2_GBPres)
        + ( (1.0 - h_star) * H_star ) * (r1_GBPmut + r2_BGPres)
        + ( (1.0 - h_star) * (1.0 - H_star) ) * (r1_BBPmut + r2_BBPres)
        - 2 * H_star;
    assert( std::abs(H_dot) < 1.0e-6 );
    return num / den;
  }
};

#endif