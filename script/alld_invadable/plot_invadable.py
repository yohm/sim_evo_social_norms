# %%
import numpy as np
import matplotlib.pyplot as plt
import os,sys

# %%
script_path = os.path.dirname(os.path.abspath(__file__))
if 'ipykernel' in sys.modules:
    # we are in a jupyter notebook
    input_file = os.path.join(script_path, '../../script_grouped/fix_prob_results/third_order_mu0.01/fixation_probs_15.dat')
else:
    # we are in a python script
    # input file is specified as a command line argument
    if len(sys.argv) < 2:
        print("Usage: python alld_invadable.py <input_file>")
        sys.exit(1)
    input_file = sys.argv[1]

input_file
# %%
# Open the file for reading

norm_ids = []
self_coop_level = []
p_fix = []


with open(input_file, 'r') as file:
    # Iterate through each line in the file
    for line in file:
        if line.startswith('#'):
            continue
        # Split the line into individual values
        values = line.split()

        # Parse the first two values and append them to the respective lists
        norm_ids.append(int(values[0]))
        self_coop_level.append(float(values[1]))

        # Parse the remaining values and append them to the data_matrix list
        p_fix.append([float(val) for val in values[2:]])

# %%
# convert lists into numpy arrays
norm_ids = np.array(norm_ids)
self_coop_level = np.array(self_coop_level)
p_fix = np.array(p_fix)

norm_ids, self_coop_level, p_fix
# %%
p_fix[0,:]
# %%
def get_invadability(norm_id):
    idx = np.where(norm_ids == norm_id)[0][0]
    invadability = p_fix[idx,:]
    # sort the invadability in descending order
    # and return the sorted invadability, norm_ids and self_coop_level
    sorted_indices = np.argsort(invadability)[::-1]
    invadability_sorted = invadability[sorted_indices]
    norm_ids_sorted = norm_ids[sorted_indices]
    self_coop_level_sorted = self_coop_level[sorted_indices]
    return invadability_sorted, norm_ids_sorted, self_coop_level_sorted

# %%
alld_invadability_sorted, norm_ids_sorted, self_coop_level_sorted = get_invadability(64704)
alld_invadability_sorted, norm_ids_sorted, self_coop_level_sorted
# %%
from matplotlib.cm import ScalarMappable

# plot a bar chart of the fixation probabilities
plt.clf()
fig, ax = plt.subplots(figsize=(6, 4))
colormap = plt.get_cmap('RdYlBu')
clormap = colormap
sm = ScalarMappable(cmap=colormap)
cbar = plt.colorbar(sm, ax=ax)
cbar.set_label('self-cooperation level\nof mutants', fontsize=12)

ax.set_xlim(-1, 81)
ax.set_ylim(0, 0.11)
for i in range(81):
    c = colormap(self_coop_level_sorted[i], alpha=0.9)
    ax.bar(i, alld_invadability_sorted[i], color=c, width=1, edgecolor='#222222', linewidth=0.3)
# plot y = 0.02
ax.plot([-1, 81], [0.02, 0.02], color='#666666', linestyle='--', linewidth=1.4)
ax.set_xlabel('rank', fontsize=16)
ax.set_ylabel('fixation probability\nagainst AllD residents', fontsize=16)
# %%
fig.savefig('alld_invadability.pdf', bbox_inches='tight', pad_inches=0.05)
# %%

l1_invadability_sorted, norm_ids_sorted, self_coop_level_sorted = get_invadability(765131)
l1_invadability_sorted, norm_ids_sorted, self_coop_level_sorted
# %%
import matplotlib.colors as mcolors

plt.clf()
fig, ax = plt.subplots(figsize=(6, 4))
colormap = plt.get_cmap('RdYlBu')
colors_with_alpha = colormap(np.arange(colormap.N))
colors_with_alpha[:, -1] = 0.9
custom_cmap = mcolors.ListedColormap(colors_with_alpha)

#clormap = colormap
sm = ScalarMappable(cmap=custom_cmap)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax)
cbar.set_label('self-cooperation level\nof mutants', fontsize=12)

xmax = 201
ax.set_xlim(-1, xmax)
ax.set_ylim(0, 0.059)
for i in range(xmax):
    c = custom_cmap(self_coop_level_sorted[i])
    ax.bar(i, l1_invadability_sorted[i], color=c, width=1.03) #, edgecolor='#222222', linewidth=0.3)
# plot y = 0.02
ax.plot([-1, xmax], [0.02, 0.02], color='#666666', linestyle='--', linewidth=1.4)
ax.set_xlabel('rank', fontsize=16)
ax.set_ylabel('fixation probability\nagainst L1 residents', fontsize=16)

# %%
fig.savefig('l1_invadability.pdf', bbox_inches='tight', pad_inches=0.05)
# %%
