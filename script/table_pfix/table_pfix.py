# %%
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import subprocess
import os

# %%
script_path = os.path.dirname(os.path.abspath(__file__))

# %%
datpath = os.path.join(script_path, '../fix_prob_results/third_order_mu0.01/fixation_probs_15.dat')
dat = np.loadtxt(datpath)
dat

# %%
norm_names = [
  [765131, "L1"],
  [634059, "L2"],
  [769226, "L3"],
  [761034, "L4"],
  [638154, "L5"],
  [629962, "L6"],
  [859333, "L7"], # or 756938
  [892101, "L8"], # or 625866
  [1047759, "ALLC"],
  [1047744, "ALLD"]
]

# %%
# find index of named norms
indices = []
for nid, name in norm_names:
  idx = np.where(dat[:,0] == nid)
  print(f'{name} {idx}')
  indices.append(idx[0][0])
indices
# %%
for i, ii in enumerate(indices):
  row = []
  for j, jj in enumerate(indices):
    pfix = dat[ii, jj+2]
    if np.abs(pfix-0.02) < 0.01:
      row.append(f"\\underline{{{pfix:.2f}}}") 
    elif pfix > 0.02:
      row.append(f"\\textbf{{{pfix:.2f}}}") 
    else:
      row.append(f"{pfix:.2f}")
  print(' & '.join(row) + ' \\\\')


# %%
