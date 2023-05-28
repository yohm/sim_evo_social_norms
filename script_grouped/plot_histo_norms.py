#%%
import matplotlib.pyplot as plt
import numpy as np
# %%
# Load data
# read "histo_norms.dat" and store it in a variable called "data"
open_file = open('histo_norms.dat', 'r')
fracs = []
c_levels = []
names = []
for line in open_file:
  # split line into list
  line = line.split()
  frac = float(line[2])
  if frac < 0.05:
    fracs.append( 1.0 - sum(fracs) )
    c_levels.append(0)
    names.append('Other')
    break
  fracs.append(frac)
  c_levels.append(float(line[3]))
  name = line[4] + '\n' + line[1]
  names.append(name)

fracs, c_levels, names
# %%
# Plot data
plt.clf()
fig, ax = plt.subplots()

colors = []
cmap = plt.get_cmap('coolwarm')
colors = [cmap(c_level) for c_level in c_levels]

# add color scale legend
plt.colorbar(plt.cm.ScalarMappable(cmap=cmap), label='Cooperation Level', ax=ax)

wedges, texts = plt.pie(fracs, labels=names, colors=colors, counterclock=False, startangle=90)
for w in wedges:
  w.set_linewidth(0.3)
  w.set_edgecolor('black')

plt.savefig('histo_norms.pdf')


# %%
