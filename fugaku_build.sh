#!/bin/bash -ex

# Fujitsu compiler with and without clang
mpiFCCpx -stdlib=libc++ -Nclang -Kfast -Kopenmp -I$HOME/data/sandbox/json/include -I$HOME/data/sandbox/eigen-3.3.7 -Iicecream-cpp -Icaravan-lib -o main_priv_sim.out main_priv_sim.cpp


# gcc-FJMPI cross compile
# . /vol0004/apps/oss/gcc-arm-11.2.1/setup-env.sh
# g++ -static -std=c++14 -O3 -fopenmp -I$HOME/data/sandbox/json/include -I$HOME/data/sandbox/eigen-3.3.7 -Iicecream-cpp -Icaravan-lib -c -o main_priv_sim.o main_priv_sim.cpp
# mpiFCCpx -o main_priv_sim.out main_priv_sim.o

