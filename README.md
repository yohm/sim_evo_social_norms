# Prerequisites

- cmake
- OpenMP
- MPI
- Eigen3
- pagmo

Install these prerequisites by homebrwe if you are on macOS.

```bash
brew install cmake
brew install libomp
brew install open-mpi
brew install eigen
brew install pagmo
```

# Build

Clone the repository with submodules.

```bash
git clone --recursive git@github.com:yohm/sim_indirect_recip_stochastic.git
cd sim_indirect_recip_stochastic
```

Make a build directory and run cmake.

```bash
mkdir build
cd build
cmake ..
cmake --build .
```

