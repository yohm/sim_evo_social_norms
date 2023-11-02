#include <iostream>
#include <fstream>
#include <chrono>
#include <regex>
#include <icecream.hpp>
#include "GroupedEvo.hpp"

#include <gtest/gtest.h>


TEST(GroupedEvo, InterGroupImitationProb) {
  double p_imit = GroupedEvo::InterGroupImitationProb(1.0, 0.0, 5.0, 1.0);
  double expected = 1.0 / (1.0 + std::exp(1.0 * (4.0 - 0.0)));
  EXPECT_NEAR(p_imit, expected, 0.002);
}

