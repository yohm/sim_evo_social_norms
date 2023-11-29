# %%
import numpy as np
import matplotlib.pyplot as plt
import os,subprocess
import json,pickle

# %%
script_path = os.path.dirname(os.path.abspath(__file__))
exe_path = os.path.join(script_path, '..', '..', 'cmake-build-release', 'main_well_mixed_evo')
input_dir_path = os.path.join(script_path, '..', 'fix_prob_results')

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
numpy_pc3 = None
fname = 'third_order_eq_pc.dat'
if os.path.exists(fname):
  with open(fname, 'rb') as f:
    numpy_pc3 = np.loadtxt(f)
numpy_pc3

# %%
if numpy_pc3 is None:
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

  numpy_pc3 = np.array(pc_all)

# %%
np.savetxt(fname, numpy_pc3)
numpy_pc3





# %%
plt.clf()
fig,ax = plt.subplots(1,1,figsize=(6,4))

color_map = plt.get_cmap('viridis')
for mu_i,mu in enumerate(mu_list):
  ax.plot(benefit_list, numpy_pc3[mu_i,:], label=f'$\mu_a={mu}$', marker='o', color=color_map(mu_i/len(mu_list)))
ax.set_xlim([1.2,5.2])
ax.set_xticks([2.0, 3.0, 4.0, 5.0])
ax.set_yticks([0.0,0.2,0.4,0.6,0.8,1.0])
ax.set_xticklabels([2, 3, 4, 5])
ax.set_yticklabels([0.0,0.2,0.4,0.6,0.8,1.0])
ax.set_ylim([0.0,1.0])
ax.legend(loc='upper right', fontsize=12)
ax.set_xlabel('benefit', fontsize=18)
ax.set_ylabel('cooperation level', fontsize=18)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# %%
fig.savefig('third_order_evo.pdf', bbox_inches='tight')

# %%
# second-order social norms
numpy_pc2 = None
fname = 'second_order_eq_pc.dat'
if os.path.exists(fname):
  with open(fname, 'rb') as f:
    numpy_pc2 = np.loadtxt(f)
numpy_pc2

# %%
if numpy_pc2 is None:
  pc_all = []
  for mu in mu_list:
    input_dir = os.path.join(input_dir_path, f'second_order_mu{mu}')
    paths = [os.path.join(input_dir, f'fixation_probs_{s}.msgpack') for s in range(0, 8)]

    pc_list = []
    for p in paths:
      out = calc_stationary(p)
      pc = equilibrium_coop_level(out)
      pc_list.append(pc)
    pc_all.append(pc_list)

  numpy_pc2 = np.array(pc_all)

# %%
np.savetxt(fname, numpy_pc2)
numpy_pc2

# %%
plt.clf()
fig,ax = plt.subplots(1,1,figsize=(6,4))

color_map = plt.get_cmap('viridis')
for mu_i,mu in enumerate(mu_list):
  ax.plot(benefit_list, numpy_pc2[mu_i,:], label=f'$\mu_a={mu}$', marker='o', color=color_map(mu_i/len(mu_list)))
ax.set_xlim([1.2,5.2])
ax.set_xticks([2.0, 3.0, 4.0, 5.0])
ax.set_yticks([0.0,0.2,0.4,0.6,0.8,1.0])
ax.set_xticklabels([2, 3, 4, 5])
ax.set_yticklabels([0.0,0.2,0.4,0.6,0.8,1.0])
ax.set_ylim([0.0,1.0])
ax.legend(loc='upper right', fontsize=12)
ax.set_xlabel('benefit', fontsize=18)
ax.set_ylabel('cooperation level', fontsize=18)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)


# %%
fig.savefig('second_order_evo.pdf', bbox_inches='tight')

# %%
dat = {}
if os.path.exists('three_species.pickle'):
  with open('three_species.pickle', 'rb') as f:
    dat = pickle.load(f)
dat
# %%
