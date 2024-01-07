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
git clone --recursive git@github.com:yohm/sim_evo_social_norms.git
cd sim_evo_social_norms
```

Make a build directory and run cmake.

```bash
mkdir cmake-build-debug
cd cmake-build-debug
cmake ..
cmake --build .
```

If you want to build in release mode, run cmake with `-DCMAKE_BUILD_TYPE=Release` option.

```bash
mkdir cmake-build-release
cd cmake-build-release
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
```

## Executables

- `main_calc_fixation_probs`: Calculate all the fixation probabilities between a set of social norms.
- `main_well_mixed_evo`: Calculate the stationary distribution of evolutionary process in well-mixed populations.
- `main_grouped_evo_ODE`: Calculate the stationary distribution of evolutionary process in group-structured populations using ODE in the M-infinity limit.
- `main_grouped_evo_MC`: Calculate the stationary distribution of evolutionary process in group-structured populations using Monte Carlo simulation with finite M.
- `main_fix_probs_param_dep`: Calculate the fixation probability between two species for different sigma & benefit.
- `inspect_Norm`, `inspect_PublicRepGame`, `inspect_PrivRepGame`, `inspect_EvolPrivRepGame`, `inspect_GroupedEvo`: Print the details of an instance of each class.
- `msgpack2json`: Convert a msgpack file to a json file.
- `extract_2nd_order_msgpack`: Extract the fixation probabilities of the second order norms from an output msgpack file of `main_calc_fixation_probs`.

## Tests

Unit tests are prepared. The executables that starts with `test_` are the unit tests. Run these tests using `ctest` command.

```bash
cd cmake-build-debug
ctest
```

## Scripts

- `fugaku_build_calc_fixation_probs.sh`: Build `main_calc_fixation_probs` on Fugaku.
- `fugaku_job`
  - `fugaku_job.sh`: A sample job script for Fugaku.
  - `params.json`: A sample parameter file for Fugaku.
- `script`
  - `evo_three_species`: scripts for the three-species evolutionary simulations.
  - `evo_well_mixed`: scripts for the evolutionary simulations in well-mixed populations.
  - `evo_grouped_ODE`: scripts for the evolutionary simulations in grouped-structured populations using ODE
  - `evo_grouped_MC`: scripts for the evolutionary simulations in grouped-structured populations using Monte Carlo simulation
