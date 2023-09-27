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

# how to run

## `main_calc_fixation_probs`

1. Run `main_calc_fixation_probs` to calculate the fixation probabilities between social norms.
  - Since its calculation is heavy, it must be run on a supercomputer.
  - To build on Fugaku, execute `./fugaku_build_calc_fixation_probs.sh`.
2. Submit a job
  - Use scripts in `fugaku_job/` directory.
  - `fugaku_job.sh` is a sample job script. The input parameters are specified by `params.json`. When you run a job, copy these files into another directory and change parameters as you like. Don't forget to correctly specify the path to the executable in the job script if you change the working directory.
3. After the job is complete, you will find files `fixation_probs_%d.dat` and `fixation_probs_%d.msgpack`.
  - These output files contain the fixation probabilities for all possible combinations of the strategies.
  - Both `dat` and `msgpack` files contain almost the same information. `dat` file is a text file that is convenient for debugging and a quick inspection. `msgpack` file is a binary file including more detailed information.

## `main_well_mixed_evo`

Calculate the stationary distribution of evolutionary process in well-mixed populations.
Using the output of `main_calc_fixation_probs`, principal eigenvector of the Markov process is calculated via the power iteration method.

1. Run with the msgpack file as an argument: `main_well_mixed_evo <msgpack file>`
2. The abundance of each strategy in equilibrium state is printed.



