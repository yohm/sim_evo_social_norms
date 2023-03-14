cmake_minimum_required(VERSION 3.7)
project(sim_indirect_recip_stochastic)

set(CMAKE_CXX_STANDARD 17)

#set(CMAKE_VERBOSE_MAKEFILE 1)

find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(MPI REQUIRED)
find_package(pagmo REQUIRED)
find_package(nlohmann_json REQUIRED)
include_directories(${CMAKE_SOURCE_DIR}/icecream-cpp)

set(SOURCE_FILES Norm.hpp PublicRepGame.hpp)

add_executable(main_public_rep.out main_public_rep.cpp Norm.hpp PublicRepGame.hpp)

add_executable(main_priv_sim.out main_priv_sim.cpp Norm.hpp PrivRepGame.hpp)
target_link_libraries(main_priv_sim.out Eigen3::Eigen MPI::MPI_CXX nlohmann_json::nlohmann_json)
target_include_directories(main_priv_sim.out PRIVATE ${CMAKE_SOURCE_DIR}/caravan-lib)

add_executable(main_optimize_norm main_optimize_norm.cpp Norm.hpp PrivRepGame.hpp)
target_link_libraries(main_optimize_norm OpenMP::OpenMP_CXX Eigen3::Eigen MPI::MPI_CXX Pagmo::pagmo)

add_executable(test_Norm test_Norm.cpp Norm.hpp)

add_executable(test_PublicRepGame test_PublicRepGame.cpp PublicRepGame.hpp Norm.hpp)

add_executable(test_PrivRepGame test_PrivRepGame.cpp PrivRepGame.hpp Norm.hpp)
target_link_libraries(test_PrivRepGame OpenMP::OpenMP_CXX Eigen3::Eigen)

add_executable(minimal_benchmark minimal_benchmark.cpp Norm.hpp PrivRepGame.hpp)
