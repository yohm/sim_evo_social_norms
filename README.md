# Prerequisites

- cmake
- OpenMP
- MPI
- Eigen3
- nlohmann-json

Install these prerequisites by homebrew if you are on macOS.

```bash
brew install cmake
brew install libomp
brew install open-mpi
brew install eigen
brew install nlohmann-json
```

# Build

Clone the repository with submodules.

```bash
git clone --recursive git@github.com:yohm/sim_indirect_stochastic_privrep.git
cd sim_indirect_recip_stochastic
```

Make a build directory and run cmake.

```bash
mkdir build
cd build
cmake ..
cmake --build .
```

## Executables

- `main_calc_fixation_probs`: Calculate all the fixation probabilities between a set of social norms.
- `main_fix_probs_param_dep`: Calculate the fixation probability between two species for different sigma & benefit.
- `main_grouped_evo`: evolutionary simulation of group-structured population. Fixation probabilities must be calculated in advance by `main_calc_fixation_probs`.

## Tests

Unit tests are prepared. The executables that starts with `test_` are the unit tests. Run these like

```bash
./test_Norm
```

