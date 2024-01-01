# %%
import os
import numpy as np
import matplotlib.pyplot as plt

# %%
# second-order strategies, mut_r = 0.01
plt.clf()
fig,ax = plt.subplots(1,1,figsize=(6,4))

benefit_list = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
mu_list = np.array([0.01, 0.02, 0.05, 0.1])
pc = np.array([
  [0.250589, 0.32414, 0.444474, 0.728104, 0.839369, 0.867824, 0.854194, 0.804408],
  [0.169904, 0.250466, 0.305859, 0.554437, 0.703768, 0.720381, 0.652992, 0.536492],
  [0.000400532, 0.11022, 0.171466, 0.189772, 0.242828, 0.268335, 0.222071, 0.155533],
  [2.75427e-05, 0.067184, 0.0536037, 0.090105, 0.0784443, 0.0653663, 0.0553425, 0.0473016]
])

color_map = plt.get_cmap('viridis')
for mu_i,mu in enumerate(mu_list):
  ax.plot(benefit_list, pc[mu_i,:], label=f'$\mu_a={mu}$', marker='o', color=color_map(mu_i/len(mu_list)))
  ax.set_xlim([1.2,5.2])
  ax.set_xticks([2.0, 3.0, 4.0, 5.0])
  ax.set_yticks([0.0,0.2,0.4,0.6,0.8,1.0])
  ax.set_xticklabels([2.0, 3.0, 4.0, 5.0], fontsize=16)
  ax.set_yticklabels([0.0,0.2,0.4,0.6,0.8,1.0], fontsize=16)
  ax.set_ylim([0.0,1.0])
  ax.legend(loc='upper left', fontsize=18)
  ax.set_xlabel('benefit', fontsize=24)
  ax.set_ylabel('cooperation level', fontsize=24)

# %%
fig.savefig("grouped_second_order_pc_r0.01.pdf", bbox_inches="tight", pad_inches=0.3)
# %%
script_path = os.path.dirname(os.path.abspath(__file__))
input_dir_path = os.path.join(script_path, '..', 'fix_prob_results')

benefit_list = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
mu_list = np.array([0.01, 0.02, 0.05, 0.1])

pc_all2 = []
for mu in mu_list:
  pc_benefit = []
  for i,benefit in enumerate(benefit_list):
    pcs = []
    for s in range(5):
      path = os.path.join(input_dir_path, f'second_order_mu{mu}_batch', f'seed{s}', f"grouped_timeseries_{i}.dat")
      dat = np.loadtxt(path)
      pc = dat[-1,1]
      pcs.append(pc)
    print(pcs)
    pc_benefit.append([benefit, np.mean(pcs), np.std(pcs)/np.sqrt(len(pcs))])
  pc_all2.append(pc_benefit)
pc_all2 = np.array(pc_all2)
pc_all2, pc_all2.shape

# %%
# second-order strategies, mut_r = 0.05
plt.clf()
fig,ax = plt.subplots(1,1,figsize=(6,4))

color_map = plt.get_cmap('viridis')
for mu_i,mu in enumerate(mu_list):
  #ax.plot(pc_all2[mu_i,:,0], pc_all2[mu_i,:,1], label=f'$\mu_a={mu}$', marker='o', color=color_map(mu_i/len(mu_list)))
  ax.errorbar(pc_all2[mu_i,:,0], pc_all2[mu_i,:,1], yerr=pc_all2[mu_i,:,2], label=f'$\mu_a={mu}$', marker='o', color=color_map(mu_i/len(mu_list)))
ax.set_xlim([1.2,5.2])
ax.set_xticks([2.0, 3.0, 4.0, 5.0])
ax.set_yticks([0.0,0.2,0.4,0.6,0.8,1.0])
ax.set_xticklabels([2, 3, 4, 5])
ax.set_yticklabels([0.0,0.2,0.4,0.6,0.8,1.0])
ax.set_ylim([0.0,1.0])
ax.legend(loc='upper left', fontsize=12)
ax.set_xlabel('benefit', fontsize=18)
ax.set_ylabel('cooperation level', fontsize=18)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
# %%
fig.savefig("grouped_second_order_pc_r0.05.pdf", bbox_inches="tight")

# %%
# third-order strategies, mut_r = 0.05
pc_all3 = []
for mu in mu_list:
  pc_benefit = []
  for i,benefit in enumerate(benefit_list):
    pcs = []
    for s in range(5):
      path = os.path.join(input_dir_path, f'third_order_mu{mu}_batch', f'seed{s}', f"grouped_timeseries_{i}.dat")
      dat = np.loadtxt(path)
      pc = dat[-1,1]
      pcs.append(pc)
    print(pcs)
    pc_benefit.append([benefit, np.mean(pcs), np.std(pcs)/np.sqrt(len(pcs))])
  pc_all3.append(pc_benefit)
pc_all3 = np.array(pc_all3)
pc_all3, pc_all3.shape

# %%
plt.clf()
fig,ax = plt.subplots(1,1,figsize=(6,4))

color_map = plt.get_cmap('viridis')
for mu_i,mu in enumerate(mu_list):
  #ax.plot(pc_all2[mu_i,:,0], pc_all2[mu_i,:,1], label=f'$\mu_a={mu}$', marker='o', color=color_map(mu_i/len(mu_list)))
  ax.errorbar(pc_all3[mu_i,:,0], pc_all3[mu_i,:,1], yerr=pc_all3[mu_i,:,2], label=f'$\mu_a={mu}$', marker='o', color=color_map(mu_i/len(mu_list)))
ax.set_xlim([1.2,5.2])
ax.set_xticks([2.0, 3.0, 4.0, 5.0])
ax.set_yticks([0.0,0.2,0.4,0.6,0.8,1.0])
ax.set_xticklabels([2, 3, 4, 5])
ax.set_yticklabels([0.0,0.2,0.4,0.6,0.8,1.0])
ax.set_ylim([0.0,1.0])
ax.legend(loc='lower right', fontsize=12)
ax.set_xlabel('benefit', fontsize=18)
ax.set_ylabel('cooperation level', fontsize=18)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
# %%
# %%
fig.savefig("grouped_third_order_pc_r0.05.pdf", bbox_inches="tight")

# %%
# plot N-dependency
pc_all_n = []
n_list = [30, 50, 70, 100]
for n in n_list:
  pc_benefit = []
  for i,benefit in enumerate(benefit_list):
    path = os.path.join(input_dir_path, f'third_order_N{n}_mu0.02', f"well_mixed_evo_{i}.dat")
    dat = np.loadtxt(path)
    pc = dat[-1,1]
    pc_benefit.append([benefit, pc])
