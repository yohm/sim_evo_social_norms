# %%
import matplotlib.pyplot as plt
import numpy as np

# %%
open_file = open('histo_norms.dat', 'r')
dat = np.loadtxt('histo_norms.dat')

names = []
fractions = []
c_levels = []
current_sum = 0.0
other_frac = 0.0
other_c_level = 0.0

for line in open_file:
  line = line.split()
  name = int(line[0])
  frac = float(line[1])
  c_level = float(line[2])

  if frac > 0.01 and current_sum < 0.8:
    names.append(name)
    fractions.append(frac)
    c_levels.append(c_level)
  else:
    other_frac += frac
    other_c_level += frac * c_level

other_c_level /= other_frac
names.append('others')
fractions.append(other_frac)
c_levels.append(other_c_level)

# %%
norm_names = {
  765131: "L1",
  634059: "L2",
  769226: "L3",
  761034: "L4",
  638154: "L5",
  629962: "L6",
  859333: "L7",  # or 756938
  625866: "L8",
  765130: "L1 BBD",
  765129: "L1 BGD",
  634058: "L2 BBD",
  769227: "L3 BBC"
}

for i, name in enumerate(names):
  if name in norm_names:
    names[i] = norm_names[name]
names
# %%
fig, ax = plt.subplots(figsize=(6,4))
colors = []
cmap = plt.get_cmap('coolwarm')
colors = [cmap(c_level) for c_level in c_levels]

print(names, fractions, c_levels)

# add color scale legend
plt.colorbar(plt.cm.ScalarMappable(cmap=cmap), label='Cooperation Level', ax=ax)

wedges, texts = plt.pie(fractions, labels=names, colors=colors, counterclock=False, startangle=90)
for w in wedges:
  w.set_linewidth(0.3)
  w.set_edgecolor('black')

fig.savefig('histo_norms.pdf', bbox_inches='tight')

# %%
