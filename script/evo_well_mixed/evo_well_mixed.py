# %%
import numpy as np
import matplotlib.pyplot as plt
import os,subprocess
import json,pickle

# %%
script_path = os.path.dirname(os.path.abspath(__file__))
exe_path = os.path.join(script_path, '..', '..', 'cmake-build-release', 'main_well_mixed_evo')
input_dir_path = os.path.join(script_path, 'fix_prob_results')

# %%
def calc_stationary(msgpack_path):
  str = subprocess.check_output([exe_path, msgpack_path], universal_newlines=True)
  lines = str.strip().split('\n')
  data = [line.split() for line in lines]
  a = np.array(data, dtype=float)
  return a

# %%
out = calc_stationary(os.path.join(input_dir_path, 'third_order_mu0.1/fixation_probs_23.msgpack'))
out

# %%
def equilibrium_coop_level(dat):
  return np.dot( dat[:,1], dat[:,2] )

# %%
pc = equilibrium_coop_level(out)
pc


# %%
benefit_list = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
mu_list = [0.01, 0.02, 0.05, 0.1]

input_dirs = [os.path.join(input_dir_path, f'third_order_mu{mu}') for mu in mu_list]

input_paths = []
for s in range(8, 16):
  msgpack_path = os.path.join(input_dir_path, f'third_order_mu0.1/fixation_probs_{s}.msgpack')
  input_paths.append(msgpack_path)
input_paths


# %%%%%%%%%%%%%%%%%%%

# %%
plt.clf()
resident = 'L1'
color_map = plt.get_cmap('viridis')
plt.xlabel('benefit')
plt.ylabel('cooperation level')
plt.xlim([1.2,5.2])
plt.ylim([0.0,1.0])
for mu_i,mu in enumerate(mu_list):
  key = (N,mu,resident)
  plt.plot(benefit_list, dat[key], label=f'$\mu_a={mu}$', marker='o', color=color_map(mu_i/len(mu_list)))
plt.legend(loc='upper left')


# %%
fig,axs = plt.subplots(2, 4, figsize=(24,12), sharex=True, sharey=True)
plt.subplots_adjust(wspace=0.02, hspace=0.2)
custom_x_ticks = [2.0, 3.0, 4.0, 5.0]

for i, (ax, resident) in enumerate(zip(axs.flat, leading_eight)):
  ax.set_title(resident, fontsize=26)
  for mu_i,mu in enumerate(mu_list):
    key = (N,mu,resident)
    ax.plot(benefit_list, dat[key], label=f'$\mu_a={mu}$', marker='o', color=color_map(mu_i/len(mu_list)))
    ax.set_xlim([1.2,5.2])
    ax.set_xticks(custom_x_ticks)
    ax.set_xticklabels(custom_x_ticks, fontsize=16)
    ax.set_yticklabels([0.0,0.2,0.4,0.6,0.8,1.0], fontsize=16)
    ax.set_ylim([0.0,1.0])
    if i == 3:
      ax.legend(loc='upper right', fontsize=20)
    if i >= 4:
      ax.set_xlabel('benefit', fontsize=24)
    if i % 4 == 0:
      ax.set_ylabel('cooperation level', fontsize=24)

# %%
fig.savefig('three_species.pdf', bbox_inches='tight', pad_inches=0.3)
# %%

fig,axs = plt.subplots(4, 4, figsize=(24,24), sharex=True, sharey=True)
plt.subplots_adjust(wspace=0.02, hspace=0.2)
custom_x_ticks = [2.0, 3.0, 4.0, 5.0]

for i, (ax, resident) in enumerate(zip(axs.flat, secondary_sixteen)):
  ax.set_title(resident, fontsize=26)
  for mu_i,mu in enumerate(mu_list):
    key = (N,mu,resident)
    ax.plot(benefit_list, dat[key], label=f'$\mu_a={mu}$', marker='o', color=color_map(mu_i/len(mu_list)))
    ax.set_xlim([1.2,5.2])
    ax.set_xticks(custom_x_ticks)
    ax.set_xticklabels(custom_x_ticks, fontsize=16)
    ax.set_yticklabels([0.0,0.2,0.4,0.6,0.8,1.0], fontsize=16)
    ax.set_ylim([0.0,1.0])
    if i == 3:
      ax.legend(loc='upper right', fontsize=20)
    if i >= 12:
      ax.set_xlabel('benefit', fontsize=24)
    if i % 4 == 0:
      ax.set_ylabel('cooperation level', fontsize=24)
# %%
fig.savefig('three_species_S16.pdf', bbox_inches='tight', pad_inches=0.3)

# %%
