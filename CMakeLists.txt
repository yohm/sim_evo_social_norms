cmake_minimum_required(VERSION 3.7)
project(sim_indirect_stochastic_privrep)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")

#set(CMAKE_VERBOSE_MAKEFILE 1)

find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(MPI REQUIRED)
find_package(nlohmann_json REQUIRED)
include_directories(${CMAKE_SOURCE_DIR}/icecream-cpp)

# See http://google.github.io/googletest/quickstart-cmake.html
include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(main_priv_sim.out main_priv_sim.cpp Norm.hpp PrivRepGame.hpp)
target_link_libraries(main_priv_sim.out Eigen3::Eigen OpenMP::OpenMP_CXX MPI::MPI_CXX nlohmann_json::nlohmann_json)
target_include_directories(main_priv_sim.out PRIVATE ${CMAKE_SOURCE_DIR}/caravan-lib)

add_executable(main_calc_fixation_probs main_calc_fixation_probs.cpp Norm.hpp PrivRepGame.hpp Vector2d.hpp Parameters.hpp)
target_link_libraries(main_calc_fixation_probs Eigen3::Eigen OpenMP::OpenMP_CXX MPI::MPI_CXX nlohmann_json::nlohmann_json)

add_executable(main_grouped_evo main_grouped_evo.cpp Norm.hpp PrivRepGame.hpp Vector2d.hpp Parameters.hpp)
target_link_libraries(main_grouped_evo Eigen3::Eigen nlohmann_json::nlohmann_json)

add_executable(test_Norm test_Norm.cpp Norm.hpp)
target_link_libraries(test_Norm GTest::gtest_main)

add_executable(test_PublicRepGame test_PublicRepGame.cpp PublicRepGame.hpp Norm.hpp)
target_link_libraries(test_PublicRepGame GTest::gtest_main)

add_executable(test_PrivRepGame test_PrivRepGame.cpp PrivRepGame.hpp Norm.hpp)
target_link_libraries(test_PrivRepGame OpenMP::OpenMP_CXX Eigen3::Eigen nlohmann_json::nlohmann_json GTest::gtest_main)

add_executable(inspect_Norm inspect_Norm.cpp Norm.hpp)

add_executable(inspect_PublicRepGame inspect_PublicRepGame.cpp PublicRepGame.hpp Norm.hpp)

add_executable(inspect_PrivRepGame inspect_PrivRepGame.cpp PrivRepGame.hpp Norm.hpp)
target_link_libraries(inspect_PrivRepGame OpenMP::OpenMP_CXX Eigen3::Eigen nlohmann_json::nlohmann_json)

add_executable(inspect_PrivRepGameTimeSeries inspect_PrivRepGameTimeSeries.cpp PrivRepGame.hpp Norm.hpp)
target_link_libraries(inspect_PrivRepGameTimeSeries Eigen3::Eigen nlohmann_json::nlohmann_json)

add_executable(inspect_finite_mu inspect_finite_mu.cpp Norm.hpp PrivRepGame.hpp)
target_link_libraries(inspect_finite_mu Eigen3::Eigen OpenMP::OpenMP_CXX nlohmann_json::nlohmann_json)

