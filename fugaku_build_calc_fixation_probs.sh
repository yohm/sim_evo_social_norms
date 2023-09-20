#!/bin/bash -ex

# Fujitsu compiler with and without clang
mpiFCCpx -std=c++17 -stdlib=libc++ -Nclang -Kfast -Kopenmp -I$HOME/data/sandbox/json/include -I$HOME/data/sandbox/eigen-3.3.7 -Iicecream-cpp -o main_calc_fixation_probs.out main_calc_fixation_probs.cpp

