# %%
import numpy as np
import matplotlib.pyplot as plt
import os,subprocess

# %%
script_path = os.path.dirname(os.path.abspath(__file__))
input_dir_path = os.path.join(script_path, '..', 'fix_prob_results')


# %%
def equilibrium_coop_level(dat):
  return np.dot( dat[:,1], dat[:,2] )


# %%
benefit_list = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
mu_list = [0.01, 0.02, 0.05, 0.1]

# %%
pc_all3 = []
for mu in mu_list:
  pc_benefit = []
  for i,benefit in enumerate(benefit_list):
    pcs = []
    for s in range(5):
      path = os.path.join(input_dir_path, f'third_order_mu{mu}_batch', f'seed{s}', f"well_mixed_evo_{i}.txt")
      dat = np.loadtxt(path)
      pc = equilibrium_coop_level(dat)
      pcs.append(pc)
    pc_benefit.append([benefit, np.mean(pcs), np.std(pcs)/np.sqrt(len(pcs))])
  pc_all3.append(pc_benefit)
pc_all3 = np.array(pc_all3)
pc_all3, pc_all3.shape

# %%
plt.clf()
fig,ax = plt.subplots(1,1,figsize=(6,4))

color_map = plt.get_cmap('viridis')
for mu_i,mu in enumerate(mu_list):
  ax.errorbar(pc_all3[mu_i,:,0], pc_all3[mu_i,:,1], yerr=pc_all3[mu_i,:,2], label=f'$\epsilon={mu}$', marker='o', color=color_map(mu_i/len(mu_list)))
  #ax.plot(pc_all3[mu_i,:,0], pc_all3[mu_i,:,1], label=f'$\mu_a={mu}$', marker='o', color=color_map(mu_i/len(mu_list)))
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
pc_all2 = []
for mu in mu_list:
  pc_benefit = []
  for i,benefit in enumerate(benefit_list):
    pcs = []
    for s in range(5):
      path = os.path.join(input_dir_path, f'second_order_mu{mu}_batch', f'seed{s}', f"well_mixed_evo_{i}.txt")
      dat = np.loadtxt(path)
      pc = equilibrium_coop_level(dat)
      pcs.append(pc)
    pc_benefit.append([benefit, np.mean(pcs), np.std(pcs)/np.sqrt(len(pcs))])
  pc_all2.append(pc_benefit)
pc_all2 = np.array(pc_all2)
pc_all2, pc_all2.shape

# %%
plt.clf()
fig,ax = plt.subplots(1,1,figsize=(6,4))

color_map = plt.get_cmap('viridis')
for mu_i,mu in enumerate(mu_list):
  ax.errorbar(pc_all2[mu_i,:,0], pc_all2[mu_i,:,1], yerr=pc_all2[mu_i,:,2], label=f'$\epsilon={mu}$', marker='o', color=color_map(mu_i/len(mu_list)))
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
# N-dependency

pc_all_n = []
n_list = [30, 50, 70, 100]
for n in n_list:
  pc_benefit = []
  for i,benefit in enumerate(benefit_list):
    path = os.path.join(input_dir_path, f'third_order_N{n}_mu0.02', f"well_mixed_evo_{i}.dat")
    dat = np.loadtxt(path)
    pc = equilibrium_coop_level(dat)
    pc_benefit.append([benefit, pc])
  pc_all_n.append(pc_benefit)
pc_all_n = np.array(pc_all_n)
pc_all_n, pc_all_n.shape
# %%
plt.clf()
fig,ax = plt.subplots(1,1,figsize=(6,4))

color_map = plt.get_cmap('plasma')
for n_i,n in enumerate(n_list):
  ax.plot(pc_all_n[n_i,:,0], pc_all_n[n_i,:,1], label=f'$N={n}$', marker='o', color=color_map(n_i/len(n_list)))
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
fig.savefig('third_order_N_dep.pdf', bbox_inches='tight')
# %%
# mu_e-dependency

pc_all_n = []
mue_list = [0, 0.01, 0.02, 0.05]
for mue in mue_list:
  pc_benefit = []
  for i,benefit in enumerate(benefit_list):
    path = os.path.join(input_dir_path, f'third_order_mu0.02_mue{mue}', f"well_mixed_evo_{i}.dat")
    print(path)
    dat = np.loadtxt(path)
    pc = equilibrium_coop_level(dat)
    pc_benefit.append([benefit, pc])
  pc_all_n.append(pc_benefit)
pc_all_n = np.array(pc_all_n)
pc_all_n, pc_all_n.shape
# %%
plt.clf()
fig,ax = plt.subplots(1,1,figsize=(6,4))

color_map = plt.get_cmap('plasma')
for mue_i,mue in enumerate(mue_list):
  ax.plot(pc_all_n[mue_i,:,0], pc_all_n[mue_i,:,1], label=f'$\epsilon_I={mue}$', marker='o', color=color_map(mue_i/len(mue_list)))
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
fig.savefig('third_order_mue_dep.pdf', bbox_inches='tight')
# %%
