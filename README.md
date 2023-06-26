# Prerequisites

- cmake
- OpenMP
- MPI
- Eigen3
- pagmo
- nlohmann-json

Install these prerequisites by homebrew if you are on macOS.

```bash
brew install cmake
brew install libomp
brew install open-mpi
brew install eigen
brew install pagmo
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

