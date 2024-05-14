
# %%
import numpy as np
import matplotlib.pyplot as plt
import os,subprocess

# %%
script_path = os.path.dirname(os.path.abspath(__file__))
exe_path = os.path.join(script_path, '..', '..', 'cmake-build-release', 'main_well_mixed_evo')
input_dir_path = os.path.join(script_path, '..', 'fix_prob_results')

# %%
def calc_stationary(msgpack_path):
  str = subprocess.check_output([exe_path, msgpack_path], universal_newlines=True)
  lines = str.strip().split('\n')
  names = []
  fracs = []
  c_levels = []
  for line in lines:
    line = line.split()
    names.append(int(line[0]))
    fracs.append(float(line[1]))
    c_levels.append(float(line[2]))
  return np.array(names), np.array(fracs), np.array(c_levels)
# %%
names, fracs, c_levels = calc_stationary(os.path.join(input_dir_path, 'third_order_mu0.01/fixation_probs_15.msgpack'))
names, fracs, c_levels
# %%
sorted_indices = np.argsort(fracs)[::-1]
names_sorted = names[sorted_indices]
fracs_sorted = fracs[sorted_indices]
c_levels_sorted = c_levels[sorted_indices]
names_sorted, fracs_sorted, c_levels_sorted

# %%
norm_names = {
  765131: "L1",
  634059: "L2",
  769226: "L3",
  761034: "L4",
  638154: "L5",
  629962: "L6",
  859333: "L7",  # or 756938
  892101: "L8", # or 625866
  765130: "L1' L1 BBD",
  765129: "L1'' L1 BGD",
  634058: "L2' BBD",
  634057: "L2'' BGD",
  634050: "- L2 GGD BBD",
  769227: "L3' L3 BBC",
  761035: "L4' BBC",
  638155: "L5' BBC",
  859341: "L7' BBC",
  568523: "- L2 GBDB"
}
for i in range(31):
  name = names_sorted[i]
  if names_sorted[i] in norm_names:
    name = norm_names[names_sorted[i]]
  if (names_sorted[i] & 0b1111) == 0:
    name = "ALLD"
  print(name, fracs_sorted[i], c_levels_sorted[i])
# %%
from matplotlib.cm import ScalarMappable
import matplotlib.patches as patches

# plot a bar chart of the fixation probabilities
plt.clf()
fig, ax = plt.subplots(figsize=(6, 4))
colormap = plt.get_cmap('RdYlBu')
clormap = colormap
sm = ScalarMappable(cmap=colormap)
cbar = plt.colorbar(sm, ax=ax)
cbar.set_label('self-cooperation level\nof strategies', fontsize=12)
cbar.ax.tick_params(labelsize=12)

ax.set_xlim(0.5, 26)
ax.set_ylim(0, 0.025)
for i in range(26):
    c = colormap(c_levels_sorted[i], alpha=0.9)
    ax.bar(i+1, fracs_sorted[i], color=c, width=1, edgecolor='#222222', linewidth=0.3)
    if names_sorted[i] in norm_names:
      name = norm_names[names_sorted[i]]
      o1 = 1.1
      o2 = 0.0004
      if len(name) == 2:
        ax.text(i+o1, fracs_sorted[i]+o2, name, fontsize=13, rotation=90, ha='center', va='bottom')
      elif name == "L2 GBDB":
        # ax.text(i+o1, fracs_sorted[i]+o2, "L2 without justified punishment", fontsize=8, rotation=90, ha='center', va='bottom')
        pass
      elif len(name) > 2 and name.startswith("L"):
        display_name = name.split()[0]
        ax.text(i+o1, fracs_sorted[i]+o2, display_name, fontsize=13, rotation=90, ha='center', va='bottom')
# plot y = 0.02
ax.plot([-1, 81], [1.0/2080, 1.0/2080], color='#666666', linestyle='--', linewidth=1.4)
ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_xticks([1, 5, 10, 15, 20, 25])
ax.set_xticklabels([1, 5, 10, 15, 20, 25])
ax.set_yticks([0.0, 0.01, 0.02])
ax.set_xlabel('rank', fontsize=16)
ax.set_ylabel('fraction', fontsize=16)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax.text(22, 0.003, 'ALLD', fontsize=13, rotation=0, ha='center', va='bottom')

# %%
fig.savefig('well_mixed_fractions.pdf', bbox_inches='tight', pad_inches=0.15)

# %%
plt.clf()
fig, ax = plt.subplots(figsize=(6, 4))
colormap = plt.get_cmap('RdYlBu')
clormap = colormap
sm = ScalarMappable(cmap=colormap)
cbar = plt.colorbar(sm, ax=ax)
cbar.set_label('self-cooperation level\nof strategies', fontsize=12)
cbar.ax.tick_params(labelsize=12)

ax.set_xlim(-1, 501)
ax.set_ylim(0, 0.025)
for i in range(501):
    c = colormap(c_levels_sorted[i], alpha=0.9)
    ax.bar(i, fracs_sorted[i], color=c, width=1.0, linewidth=0.0)
# plot y = 0.02
ax.tick_params(axis='both', which='major', labelsize=16)
ax.set_xticks([1, 250, 500])
ax.set_yticks([0.0, 0.01, 0.02])
ax.plot([-1, 501], [1.0/2080, 1.0/2080], color='#666666', linestyle='--', linewidth=1.4)
ax.text(100, 0.003, 'ALLD', fontsize=13, rotation=0, ha='center', va='bottom')
ax.set_xlabel('rank', fontsize=16)
ax.set_ylabel('fraction', fontsize=16)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)


# %%
fig.savefig("well_mixed_fractions_all.pdf", bbox_inches="tight", pad_inches=0.15)

# %%
names, fracs, c_levels = calc_stationary(os.path.join(input_dir_path, 'second_order_mu0.01/fixation_probs_7.msgpack'))
names, fracs, c_levels
# %%
sorted_indices = np.argsort(fracs)[::-1]
names_sorted = names[sorted_indices]
fracs_sorted = fracs[sorted_indices]
c_levels_sorted = c_levels[sorted_indices]
names_sorted, fracs_sorted, c_levels_sorted

# %%
for i in range(36):
  name = names_sorted[i]
  if names_sorted[i] in norm_names:
    name = norm_names[names_sorted[i]]
  elif (names_sorted[i] & 0b1111) == 0:
    name = "ALLD"
  elif (names_sorted[i] & 0b1111) == 0b1111:
    name = "ALLC"
  print(name, fracs_sorted[i], c_levels_sorted[i])

# %%
plt.clf()
fig, ax = plt.subplots(figsize=(6, 4))
colormap = plt.get_cmap('RdYlBu')
clormap = colormap
sm = ScalarMappable(cmap=colormap)
cbar = plt.colorbar(sm, ax=ax)
cbar.set_label('self-cooperation level\nof strategies', fontsize=12)
cbar.ax.tick_params(labelsize=12)

ax.set_xlim(0.5, 26)
ax.set_ylim(0, 0.1)
for i in range(36):
    c = colormap(c_levels_sorted[i], alpha=0.9)
    ax.bar(i+1, fracs_sorted[i], color=c, width=1, edgecolor='#222222', linewidth=0.3)
    if names_sorted[i] in norm_names:
      name = norm_names[names_sorted[i]]
      o1 = 1.1
      o2 = 0.002
      if len(name) == 2:
        ax.text(i+o1, fracs_sorted[i]+o2, name, fontsize=12, rotation=90, ha='center', va='bottom')
      elif name == "L2 GBDB":
        # ax.text(i+o1, fracs_sorted[i]+o2, "L2 without justified punishment", fontsize=8, rotation=90, ha='center', va='bottom')
        pass
      elif len(name) > 2 and name.startswith("L"):
        name = name[:2] + "'"
        ax.text(i+o1, fracs_sorted[i]+o2, name, fontsize=12, rotation=90, ha='center', va='bottom')
# plot y = 0.02
ax.plot([-1, 81], [1.0/36, 1.0/36], color='#666666', linestyle='--', linewidth=1.4)
ax.tick_params(axis='both', which='major', labelsize=16)
ax.set_xticks([1, 5, 10, 15, 20, 25])
ax.set_xticklabels([1, 5, 10, 15, 20, 25])
ax.set_xlabel('rank', fontsize=16)
ax.set_ylabel('fraction', fontsize=16)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax.text(5.5, 0.087, 'ALLD', fontsize=12, rotation=0, ha='center', va='bottom')
# %%
fig.savefig('well_mixed_fractions_second.pdf', bbox_inches='tight', pad_inches=0.15)
# %%
