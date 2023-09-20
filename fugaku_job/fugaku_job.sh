#!/bin/sh
#PJM --rsc-list "node=1200"
#PJM --rsc-list "elapse=06:00:00"
#PJM --rsc-list "rscgrp=large"
#PJM --mpi "max-proc-per-node=12"
#PJM -x PJM_LLIO_GFSCACHE=/vol0004
#PJM -S

export OMP_NUM_THREADS=4

mpiexec -stdout-proc ./%/1000R/stdout -stderr-proc ./%/1000R/stderr ../main_calc_fixation_probs.out -j params.json --norm-set third_order
#mpiexec -stdout-proc ./%/1000R/stdout -stderr-proc ./%/1000R/stderr ../main_calc_fixation_probs.out -j params.json --norm-set dual_second_order
