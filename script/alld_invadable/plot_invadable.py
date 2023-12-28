# %%
import numpy as np
import matplotlib.pyplot as plt
import os,sys
from matplotlib.cm import ScalarMappable

# %%
script_path = os.path.dirname(os.path.abspath(__file__))
if 'ipykernel' in sys.modules:
    # we are in a jupyter notebook
    input_file = os.path.join(script_path, '../fix_prob_results/third_order_mu0.01/fixation_probs_15.dat')
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
def plot_invadability(norm_id, x_max=81, y_max=0.11, bar_linewidth=0.3):
    invadability_sorted, norm_ids_sorted, self_coop_level_sorted = get_invadability(norm_id)
    # plot a bar chart of the fixation probabilities
    plt.clf()
    fig, ax = plt.subplots(figsize=(6, 4))
    colormap = plt.get_cmap('RdYlBu')
    sm = ScalarMappable(cmap=colormap)
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label('self-cooperation level\nof mutants', fontsize=12)

    ax.set_xlim(0, x_max)
    ax.set_ylim(0, y_max)
    for i in range(x_max):
        c = colormap(self_coop_level_sorted[i], alpha=0.9)
        ax.bar(i+1, invadability_sorted[i], color=c, width=1, edgecolor='#222222', linewidth=bar_linewidth)
    # plot y = 0.02
    ax.plot([-1, x_max], [0.02, 0.02], color='#666666', linestyle='--', linewidth=1.4)
    ax.set_xlabel('rank', fontsize=18)
    ax.set_ylabel('fixation probability', fontsize=18)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    return fig, ax, norm_ids_sorted

# %%
def show_norms(norm_ids_s):
    norm_names = {
        765131: "L1",
        634059: "L2",
        769226: "L3",
        761034: "L4",
        638154: "L5",
        629962: "L6",
        859333: "L7",  # or 756938
        892101: "L8", # or 625866
        765130: "L1 BBD",
        765129: "L1 BGD",
        634058: "L2 BBD",
        634057: "L2 BGD",
        634050: "L2 GGD BBD",
        769227: "L3 BBC",
        761035: "L4 BBC",
        638155: "L5 BBC",
        859341: "L7 BBC",
        568523: "L2 GBDB",
        699594: "IS"
    }
    for i,nid in enumerate(norm_ids_s):
        if i > 201:
            break
        name = nid
        if nid in norm_names:
            name = norm_names[nid]
        if (nid & 0b1111) == 0:
            name = "ALLD"
        if (nid & 0b1111) == 0b1111:
            name = "ALLC"
        print(name, i)
# %%
fig, ax, norm_ids_s = plot_invadability(64704, x_max=81, y_max=0.11, bar_linewidth=0.3)
#ax.set_ylabel(r"fixation probability, $\rho_{\rm{ALLD} \to X}$", fontsize=16)
ax.annotate('L1', xy=(1.5, 0.1015), xytext=(15, 0.103), fontsize=12, arrowprops=dict(arrowstyle='-', color='#222222'))
ax.annotate('L2', xy=(10.5, 0.0745), xytext=(24, 0.085), fontsize=12, arrowprops=dict(arrowstyle='-', color='#222222'))
ax.annotate('L7', xy=(31.5, 0.049), xytext=(44, 0.06), fontsize=12, arrowprops=dict(arrowstyle='-', color='#222222'))
ax.annotate('IS', xy=(40.5, 0.026), xytext=(52, 0.04), fontsize=12, arrowprops=dict(arrowstyle='-', color='#222222'))
ax.text(70, 0.025, 'ALLD', fontsize=12, horizontalalignment='center', verticalalignment='center')
fig.savefig('alld_invadability.pdf', bbox_inches='tight')
show_norms(norm_ids_s)
# %%
fig, ax, norm_ids_s = plot_invadability(765131, x_max=201, y_max=0.149, bar_linewidth=0.0)
ax.annotate('L7', xy=(175.5, 0.02), xytext=(150, 0.045), fontsize=12, arrowprops=dict(arrowstyle='-', color='#222222'))
ax.annotate('L4', xy=(178.5, 0.02), xytext=(168, 0.045), fontsize=12, arrowprops=dict(arrowstyle='-', color='#222222'))
ax.annotate('L3', xy=(179.5, 0.02), xytext=(186, 0.045), fontsize=12, arrowprops=dict(arrowstyle='-', color='#222222'))
ax.text(100, 0.035, 'ALLC', fontsize=12, horizontalalignment='center', verticalalignment='center')
fig.savefig('l1_invadability.pdf', bbox_inches='tight')
show_norms(norm_ids_s)

# %%
fig, ax, norm_ids_s = plot_invadability(634059, x_max=201, y_max=0.149, bar_linewidth=0.0)
ax.annotate('L1', xy=(1.5, 0.141), xytext=(50, 0.135), fontsize=12, arrowprops=dict(arrowstyle='-', color='#222222'))
ax.annotate('L7', xy=(3.5, 0.138), xytext=(50, 0.125), fontsize=12, arrowprops=dict(arrowstyle='-', color='#222222'))
ax.annotate('L3', xy=(4.5, 0.126), xytext=(50, 0.11), fontsize=12, arrowprops=dict(arrowstyle='-', color='#222222'))
ax.annotate('L4', xy=(7.5, 0.120), xytext=(50, 0.10), fontsize=12, arrowprops=dict(arrowstyle='-', color='#222222'))
ax.annotate('L5', xy=(36.5, 0.0180), xytext=(65, 0.04), fontsize=12, arrowprops=dict(arrowstyle='-', color='#222222'))
ax.text(100, 0.007, 'ALLC', fontsize=12, horizontalalignment='center', verticalalignment='center')
fig.savefig('l2_invadability.pdf', bbox_inches='tight')
show_norms(norm_ids_s)

# %%
fig, ax, norm_ids_s = plot_invadability(769226, x_max=201, y_max=0.149, bar_linewidth=0.0)
fig.savefig('l3_invadability.pdf', bbox_inches='tight', pad_inches=0.05)
show_norms(norm_ids_s)

# %%
fig, ax, norm_ids_s = plot_invadability(761034, x_max=201, y_max=0.149, bar_linewidth=0.0)
fig.savefig('l4_invadability.pdf', bbox_inches='tight', pad_inches=0.05)
show_norms(norm_ids_s)

# %%
fig, ax, norm_ids_s = plot_invadability(638154, x_max=201, y_max=0.149, bar_linewidth=0.0)
fig.savefig('l5_invadability.pdf', bbox_inches='tight', pad_inches=0.05)
show_norms(norm_ids_s)

# %%
fig, ax, norm_ids_s = plot_invadability(629962, x_max=201, y_max=1.0, bar_linewidth=0.0)
fig.savefig('l6_invadability.pdf', bbox_inches='tight', pad_inches=0.05)
show_norms(norm_ids_s)

# %%
fig, ax, norm_ids_s = plot_invadability(859333, x_max=201, y_max=0.149, bar_linewidth=0.0)
fig.savefig('l7_invadability.pdf', bbox_inches='tight', pad_inches=0.05)
show_norms(norm_ids_s)

# %%
fig, ax, norm_ids_s = plot_invadability(892101, x_max=201, y_max=0.149, bar_linewidth=0.0)
fig.savefig('l8_invadability.pdf', bbox_inches='tight', pad_inches=0.05)
show_norms(norm_ids_s)

# %%
# %%
input_file = os.path.join(script_path, '../fix_prob_results/second_order_mu0.01/fixation_probs_7.dat')
input_file

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

norm_ids = np.array(norm_ids)
self_coop_level = np.array(self_coop_level)
p_fix = np.array(p_fix)

norm_ids, self_coop_level, p_fix
# %%
fig, ax, norm_ids_s = plot_invadability(212160, x_max=36, y_max=0.11, bar_linewidth=0.3)
#ax.set_ylabel(r"fixation probability, $\rho_{\rm{ALLD} \to X}$", fontsize=16)
ax.annotate('IS', xy=(1, 0.026), xytext=(3, 0.06), fontsize=12, arrowprops=dict(arrowstyle='-', color='#222222'))
ax.annotate('L3', xy=(15, 0.0016), xytext=(20, 0.04), fontsize=12, arrowprops=dict(arrowstyle='-', color='#222222'))
ax.annotate('L6', xy=(17, 0.0), xytext=(24, 0.03), fontsize=12, arrowprops=dict(arrowstyle='-', color='#222222'))
ax.text(8, 0.02, 'ALLD', fontsize=12, horizontalalignment='center', verticalalignment='bottom')
show_norms(norm_ids_s)
fig.savefig('alld_invadability_second.pdf', bbox_inches='tight')
# %%
