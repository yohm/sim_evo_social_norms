# %%
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
# second-order strategies, mut_r = 0.05
plt.clf()
fig,ax = plt.subplots(1,1,figsize=(6,4))

benefit_list = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
mu_list = np.array([0.01, 0.02, 0.05, 0.1])
pc2 = np.array([
  [0.197408, 0.284078, 0.268917, 0.292942, 0.427655, 0.589413, 0.689211, 0.743454],
  [0.122745, 0.224372, 0.202154, 0.201286, 0.244451, 0.331971, 0.412645, 0.453978],
  [0.00165803, 0.110404, 0.128321, 0.107491, 0.094916, 0.0855022, 0.0770903, 0.0693659],
  [0.00013587, 0.0037096, 0.0357097, 0.0533747, 0.0368371, 0.0281967, 0.0237399, 0.0217193]
])

color_map = plt.get_cmap('viridis')
for mu_i,mu in enumerate(mu_list):
  ax.plot(benefit_list, pc2[mu_i,:], label=f'$\mu_a={mu}$', marker='o', color=color_map(mu_i/len(mu_list)))
ax.set_xlim([1.2,5.2])
ax.set_xticks([2.0, 3.0, 4.0, 5.0])
ax.set_yticks([0.0,0.2,0.4,0.6,0.8,1.0])
ax.set_xticklabels([2.0, 3.0, 4.0, 5.0], fontsize=16)
ax.set_yticklabels([0.0,0.2,0.4,0.6,0.8,1.0], fontsize=16)
ax.set_ylim([0.0,1.0])
ax.legend(loc='upper left', fontsize=14)
ax.set_xlabel('benefit', fontsize=18)
ax.set_ylabel('cooperation level', fontsize=18)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
# %%
fig.savefig("grouped_second_order_pc_r0.05.pdf", bbox_inches="tight", pad_inches=0.3)

# %%
# third-order strategies, mut_r = 0.05
plt.clf()
fig,ax = plt.subplots(1,1,figsize=(6,4))

benefit_list = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
mu_list = np.array([0.01, 0.02, 0.05, 0.1])
pc3 = np.array([
  [0.957341, 0.973992, 0.974725, 0.97544, 0.97586, 0.976021, 0.975982, 0.975765],
  [0.521255, 0.951134, 0.951185, 0.950901, 0.950371, 0.949537, 0.948235, 0.945916],
  [0.740454, 0.838092, 0.828762, 0.830188, 0.830839, 0.828718, 0.806295, 0.797611],
  [4.9167e-05, 0.634693, 0.715741, 0.590904, 0.573691, 0.576475, 0.584898, 0.649576]
])

color_map = plt.get_cmap('viridis')
for mu_i,mu in enumerate(mu_list):
  ax.plot(benefit_list, pc3[mu_i,:], label=f'$\mu_a={mu}$', marker='o', color=color_map(mu_i/len(mu_list)))
ax.set_xlim([1.2,5.2])
ax.set_xticks([2.0, 3.0, 4.0, 5.0])
ax.set_yticks([0.0,0.2,0.4,0.6,0.8,1.0])
ax.set_xticklabels([2.0, 3.0, 4.0, 5.0], fontsize=16)
ax.set_yticklabels([0.0,0.2,0.4,0.6,0.8,1.0], fontsize=16)
ax.set_ylim([0.0,1.0])
ax.legend(loc='lower right', fontsize=14)
ax.set_xlabel('benefit', fontsize=18)
ax.set_ylabel('cooperation level', fontsize=18)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# %%
fig.savefig("grouped_third_order_pc_r0.05.pdf", bbox_inches="tight", pad_inches=0.3)

# %%
plt.clf()
fig, axs = plt.subplots(1,2,figsize=(12,4))

for ax,pc in zip(axs, [pc2, pc3]):
  color_map = plt.get_cmap('viridis')
  for mu_i,mu in enumerate(mu_list):
    ax.plot(benefit_list, pc[mu_i,:], label=f'$\mu_a={mu}$', marker='o', color=color_map(mu_i/len(mu_list)))
  ax.set_xlim([1.2,5.2])
  ax.set_xticks([2.0, 3.0, 4.0, 5.0])
  ax.set_yticks([0.0,0.2,0.4,0.6,0.8,1.0])
  ax.set_xticklabels([2.0, 3.0, 4.0, 5.0], fontsize=16)
  ax.set_yticklabels([0.0,0.2,0.4,0.6,0.8,1.0], fontsize=16)
  ax.set_ylim([0.0,1.0])
  ax.set_xlabel('benefit', fontsize=24)
  ax.spines['right'].set_visible(False)
  ax.spines['top'].set_visible(False)

axs[0].set_ylabel('cooperation level', fontsize=18)
axs[1].legend(loc='lower right', fontsize=14)
axs[0].set_title('second-order strategies', fontsize=18)
axs[1].set_title('third-order strategies', fontsize=18)
axs[1].set_yticklabels([])
# %%
fig.savefig("grouped_second_third_r0.05.pdf", bbox_inches="tight", pad_inches=0.3)

# %%
