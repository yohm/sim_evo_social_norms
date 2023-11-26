#!/bin/bash
set -eux

script_dir=$(cd $(dirname $0); pwd)
cd $script_dir

# mu=0.01
../../../cmake-build-release/main_grouped_evo_ODE ../../fix_prob_results/second_order_mu0.01/fixation_probs_1.msgpack 1.0 0.05 500 1.0
mv timeseries.dat result/2nd_mu0.01_b2_timeseries.dat

../../../cmake-build-release/main_grouped_evo_ODE ../../fix_prob_results/second_order_mu0.01/fixation_probs_3.msgpack 1.0 0.05 500 1.0
mv timeseries.dat result/2nd_mu0.01_b3_timeseries.dat

../../../cmake-build-release/main_grouped_evo_ODE ../../fix_prob_results/second_order_mu0.01/fixation_probs_5.msgpack 1.0 0.05 500 1.0
mv timeseries.dat result/2nd_mu0.01_b4_timeseries.dat

../../../cmake-build-release/main_grouped_evo_ODE ../../fix_prob_results/second_order_mu0.01/fixation_probs_7.msgpack 1.0 0.05 500 1.0
mv timeseries.dat result/2nd_mu0.01_b5_timeseries.dat

# mu=0.02
../../../cmake-build-release/main_grouped_evo_ODE ../../fix_prob_results/second_order_mu0.02/fixation_probs_1.msgpack 1.0 0.05 500 1.0
mv timeseries.dat result/2nd_mu0.02_b2_timeseries.dat

../../../cmake-build-release/main_grouped_evo_ODE ../../fix_prob_results/second_order_mu0.02/fixation_probs_3.msgpack 1.0 0.05 500 1.0
mv timeseries.dat result/2nd_mu0.02_b3_timeseries.dat

../../../cmake-build-release/main_grouped_evo_ODE ../../fix_prob_results/second_order_mu0.02/fixation_probs_5.msgpack 1.0 0.05 500 1.0
mv timeseries.dat result/2nd_mu0.02_b4_timeseries.dat

../../../cmake-build-release/main_grouped_evo_ODE ../../fix_prob_results/second_order_mu0.02/fixation_probs_7.msgpack 1.0 0.05 500 1.0
mv timeseries.dat result/2nd_mu0.02_b5_timeseries.dat

# mu=0.05
../../../cmake-build-release/main_grouped_evo_ODE ../../fix_prob_results/second_order_mu0.05/fixation_probs_1.msgpack 1.0 0.05 500 1.0
mv timeseries.dat result/2nd_mu0.05_b2_timeseries.dat

../../../cmake-build-release/main_grouped_evo_ODE ../../fix_prob_results/second_order_mu0.05/fixation_probs_3.msgpack 1.0 0.05 500 1.0
mv timeseries.dat result/2nd_mu0.05_b3_timeseries.dat

../../../cmake-build-release/main_grouped_evo_ODE ../../fix_prob_results/second_order_mu0.05/fixation_probs_5.msgpack 1.0 0.05 500 1.0
mv timeseries.dat result/2nd_mu0.05_b4_timeseries.dat

../../../cmake-build-release/main_grouped_evo_ODE ../../fix_prob_results/second_order_mu0.05/fixation_probs_7.msgpack 1.0 0.05 500 1.0
mv timeseries.dat result/2nd_mu0.05_b5_timeseries.dat

# mu=0.1
../../../cmake-build-release/main_grouped_evo_ODE ../../fix_prob_results/second_order_mu0.1/fixation_probs_1.msgpack 1.0 0.05 500 1.0
mv timeseries.dat result/2nd_mu0.1_b2_timeseries.dat

../../../cmake-build-release/main_grouped_evo_ODE ../../fix_prob_results/second_order_mu0.1/fixation_probs_3.msgpack 1.0 0.05 500 1.0
mv timeseries.dat result/2nd_mu0.1_b3_timeseries.dat

../../../cmake-build-release/main_grouped_evo_ODE ../../fix_prob_results/second_order_mu0.1/fixation_probs_5.msgpack 1.0 0.05 500 1.0
mv timeseries.dat result/2nd_mu0.1_b4_timeseries.dat

../../../cmake-build-release/main_grouped_evo_ODE ../../fix_prob_results/second_order_mu0.1/fixation_probs_7.msgpack 1.0 0.05 500 1.0
mv timeseries.dat result/2nd_mu0.1_b5_timeseries.dat
