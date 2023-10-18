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
#out = calc_stationary(os.path.join(input_dir_path, 'third_order_mu0.1/fixation_probs_23.msgpack'))
#out

# %%
def equilibrium_coop_level(dat):
  return np.dot( dat[:,1], dat[:,2] )

# %%
#pc = equilibrium_coop_level(out)
#pc

# %%
benefit_list = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
mu_list = [0.01, 0.02, 0.05, 0.1]


# %%
numpy_pc = None
fname = 'third_order_eq_pc.dat'
if os.path.exists(fname):
  with open(fname, 'rb') as f:
    numpy_pc = np.loadtxt(f)
numpy_pc

# %%
if numpy_pc is None:
  pc_all = []
  for mu in mu_list:
    input_dir = os.path.join(input_dir_path, f'third_order_mu{mu}')
    paths = [os.path.join(input_dir, f'fixation_probs_{s}.msgpack') for s in range(8, 16)]

    pc_list = []
    for p in paths:
      out = calc_stationary(p)
      pc = equilibrium_coop_level(out)
      pc_list.append(pc)
    pc_all.append(pc_list)

  numpy_pc = np.array(pc_all)

# %%
np.savetxt(fname, numpy_pc)
numpy_pc





# %%
plt.clf()
fig,ax = plt.subplots(1,1,figsize=(8,6))

color_map = plt.get_cmap('viridis')
for mu_i,mu in enumerate(mu_list):
  ax.plot(benefit_list, numpy_pc[mu_i,:], label=f'$\mu_a={mu}$', marker='o', color=color_map(mu_i/len(mu_list)))
  ax.set_xlim([1.2,5.2])
  ax.set_xticks([2.0, 3.0, 4.0, 5.0])
  ax.set_xticklabels([2.0, 3.0, 4.0, 5.0], fontsize=16)
  ax.set_yticklabels([0.0,0.2,0.4,0.6,0.8,1.0], fontsize=16)
  ax.set_ylim([0.0,1.0])
  ax.legend(loc='upper right', fontsize=20)
  ax.set_xlabel('benefit', fontsize=24)
  ax.set_ylabel('cooperation level', fontsize=24)

fig.show()


# %%
fig.savefig('third_order_evo.pdf', bbox_inches='tight', pad_inches=0.3)
# %%
dat = {}
if os.path.exists('three_species.pickle'):
  with open('three_species.pickle', 'rb') as f:
    dat = pickle.load(f)
dat
# %%
plt.clf()
resident = 'L1'
color_map = plt.get_cmap('viridis')
plt.xlabel('benefit')
plt.ylabel('cooperation level')
plt.xlim([1.2,5.2])
plt.ylim([0.0,1.0])
for mu_i,mu in enumerate(mu_list):
  key = (50,mu,resident)
  plt.plot(benefit_list, dat[key], label=f'$\mu_a={mu}$', marker='o', color=color_map(mu_i/len(mu_list)))
plt.legend(loc='upper left')
# %%

plt.clf()
fig,axs = plt.subplots(1,2,figsize=(16,6), sharey=True)
plt.subplots_adjust(wspace=0.05)

for mu_i,mu in enumerate(mu_list):
  key = (50,mu,resident)
  axs[0].plot(benefit_list, dat[key], label=f'$\mu_a={mu}$', marker='o', color=color_map(mu_i/len(mu_list)))
axs[0].set_xticks([2.0, 3.0, 4.0, 5.0])
axs[0].set_xticklabels([2.0, 3.0, 4.0, 5.0], fontsize=16)
axs[0].set_yticks([0.0,0.2,0.4,0.6,0.8,1.0])
axs[0].set_yticklabels([0.0,0.2,0.4,0.6,0.8,1.0], fontsize=16)
axs[0].set_ylim([0.0,1.0])
axs[0].set_xlabel('benefit', fontsize=24)
axs[0].set_ylabel('cooperation level', fontsize=24)
axs[0].set_title('L1-AllC-AllD', fontsize=24)

for mu_i,mu in enumerate(mu_list):
  axs[1].plot(benefit_list, numpy_pc[mu_i,:], label=f'$\mu_a={mu}$', marker='o', color=color_map(mu_i/len(mu_list)))
axs[1].set_xlim([1.2,5.2])
axs[1].set_xticks([2.0, 3.0, 4.0, 5.0])
axs[1].set_xticklabels([2.0, 3.0, 4.0, 5.0], fontsize=16)
axs[1].set_yticks([0.0,0.2,0.4,0.6,0.8,1.0])
axs[1].set_yticklabels([0.0,0.2,0.4,0.6,0.8,1.0], fontsize=16)
axs[1].set_ylim([-0.02,1.02])
axs[1].legend(loc='upper right', fontsize=16)
axs[1].set_xlabel('benefit', fontsize=24)
axs[1].set_title('all third-order', fontsize=24)

axs[0].text(0.04, 0.95, "(a)", transform=axs[0].transAxes, fontsize=20, va='top', ha='left')
axs[1].text(0.04, 0.95, "(b)", transform=axs[1].transAxes, fontsize=20, va='top', ha='left')

# %%
fig.savefig('L1_third_order_evo.pdf', bbox_inches='tight', pad_inches=0.3)

# %%
