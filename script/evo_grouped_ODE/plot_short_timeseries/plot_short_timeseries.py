#%%
import matplotlib.pyplot as plt
import numpy as np
# %%
# Load data
dat = np.loadtxt('3rd_mu0.02_b3_timeseries.dat')

plt.clf()
fig, ax = plt.subplots(figsize=(6,4))
ax.plot(dat[:,0], dat[:,1], label='Cooperation Level')
ax.set_xlabel('time', fontsize=18)
ax.set_ylabel('cooperation level', fontsize=18)
ax.set_ylim(-0.03,1.03)

# %%
# Load data from norms.txt
# first column is ID and the remaining columns are the names
norms = []
with open('3rd_norms.txt', 'r') as file:
  line_number = 0
  for line in file:
    parts = line.strip().split(' ')
    if len(parts) >= 2:
      norms.append((int(parts[0]), ' '.join(parts[1:])))
norms

# %%
alld_idxs = []
allc_idxs = []
idx_name = {}
for idx,(nid,name) in enumerate(norms):
  if name == 'AllD':
    alld_idxs.append(idx)
  if name == 'AllC':
    allc_idxs.append(idx)
  norm_names = {
    765131: "L1",
    634059: "L2",
    769226: "L3",
    761034: "L4",
    638154: "L5",
    629962: "L6",
    859333: "L7", # or 756938
    892101: "L8", # or 625866
    765130: r"L1 with $P(B,B)=D$",
    765129: "L1 BGD",
    634058: "L2 BBD",
    769227: "L3 BBC",
    761035: "L4 BBD",
    859341: "L7 BBD"
  }
  if nid in norm_names:
    idx_name[idx] = norm_names[nid]
alld_idxs, allc_idxs, idx_name

# %%
plt.clf()
fig, ax = plt.subplots(figsize=(6,4))
ax.plot(dat[:,0], dat[:,1], label='Cooperation Level')
for idx,name in idx_name.items():
  if idx+2 < dat.shape[1]:
    if np.max(dat[:,idx+2]) > 0.05:
      ax.plot(dat[:,0], dat[:,idx+2], label=name)
# sum of ALLD
alld = np.zeros(dat.shape[0])
for idx in alld_idxs:
  if idx+2 < dat.shape[1]:
    alld += dat[:,idx+2]
allc = np.zeros(dat.shape[0])
for idx in allc_idxs:
  if idx+2 < dat.shape[1]:
    allc += dat[:,idx+2]
if np.max(alld) > 0.05:
  ax.plot(dat[:,0], alld, label='ALLD')

ax.set_xlim(0,400)
ax.set_xticks([0,100,200,300,400])
ax.set_xticklabels(['0','100','200','300','400'])
ax.set_xlabel('time', fontsize=18)
ax.set_ylabel('fractions', fontsize=18)
ax.set_ylim(-0.005,1.005)
ax.legend()

# %%
fig.savefig('3rd_timeseries.pdf', bbox_inches='tight')
# %%
dat = np.loadtxt('2nd_mu0.02_b3_timeseries.dat')

plt.clf()
fig, ax = plt.subplots(figsize=(6,4))
ax.plot(dat[:,0], dat[:,1], label='Cooperation Level')
ax.set_xlabel('time', fontsize=18)
ax.set_ylabel('cooperation level', fontsize=18)
ax.set_ylim(-0.03,1.03)
# %%
norms = []
with open('2nd_norms.txt', 'r') as file:
  line_number = 0
  for line in file:
    parts = line.strip().split(' ')
    if len(parts) >= 2:
      norms.append((int(parts[0]), ' '.join(parts[1:])))
norms

alld_idxs = []
allc_idxs = []
idx_name = {}
for idx,(nid,name) in enumerate(norms):
  if name == 'AllD':
    alld_idxs.append(idx)
  if name == 'AllC':
    allc_idxs.append(idx)
  norm_names = {
    765131: "L1",
    634059: "L2",
    769226: "L3",
    761034: "L4",
    638154: "L5",
    629962: "L6",
    859333: "L7", # or 756938
    892101: "L8", # or 625866
    765130: r"L1 with $P(B,B)=D$",
    765129: "L1 BGD",
    634058: "L2 BBD",
    769227: "L3 BBC",
    761035: "L4 BBD",
    859341: "L7 BBD"
  }
  if nid in norm_names:
    idx_name[idx] = norm_names[nid]
alld_idxs, allc_idxs, idx_name
# %%
plt.clf()
fig, ax = plt.subplots(figsize=(6,4))
ax.plot(dat[:,0], dat[:,1], label='Cooperation Level')
for idx,name in idx_name.items():
  if idx+2 < dat.shape[1]:
    if np.max(dat[:,idx+2]) > 0.05:
      ax.plot(dat[:,0], dat[:,idx+2], label=name, color=cmap(4))
# sum of ALLD
alld = np.zeros(dat.shape[0])
for idx in alld_idxs:
  if idx+2 < dat.shape[1]:
    alld += dat[:,idx+2]
allc = np.zeros(dat.shape[0])
for idx in allc_idxs:
  if idx+2 < dat.shape[1]:
    allc += dat[:,idx+2]

cmap = plt.get_cmap("tab10")
if np.max(alld) > 0.05:
  # using 3rd color in the palette
  ax.plot(dat[:,0], alld, label='ALLD', color=cmap(3))

ax.set_xlim(0,400)
ax.set_xticks([0,100,200,300,400])
ax.set_xticklabels(['0','100','200','300','400'])
ax.set_xlabel('time', fontsize=18)
ax.set_ylabel('fractions', fontsize=18)
ax.set_ylim(-0.005,1.005)
ax.legend()
# %%
fig.savefig('2nd_timeseries.pdf', bbox_inches='tight')
# %%
