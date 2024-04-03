#%%
import matplotlib.pyplot as plt
import numpy as np
# %%
# Load data
dat = np.loadtxt('result/3rd_mu0.02_b3_timeseries.dat')

plt.clf()
fig, ax = plt.subplots(figsize=(6,4))
ax.plot(dat[:,0], dat[:,1], label='Cooperation Level')
ax.set_xlabel('time', fontsize=18)
ax.set_ylabel('cooperation level', fontsize=18)
ax.set_ylim(-0.03,1.03)

# %%
# Load data from norms.txt
# first column is ID and the remaining columns are the names
norms_3rd = []
with open('3rd_norms.txt', 'r') as file:
  line_number = 0
  for line in file:
    parts = line.strip().split(' ')
    if len(parts) >= 2:
      norms_3rd.append((int(parts[0]), ' '.join(parts[1:])))
norms_3rd

# %%
alld_idxs_3rd = []
allc_idxs_3rd = []
idx_name_3rd = {}
for idx,(nid,name) in enumerate(norms_3rd):
  if name == 'AllD':
    alld_idxs_3rd.append(idx)
  if name == 'AllC':
    allc_idxs_3rd.append(idx)
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
    769227: r"L3 with $P(B,B)=C$",
    761035: "L4 BBD",
    859341: r"L7 with $P(B,B)=D$"
  }
  if nid in norm_names:
    idx_name_3rd[idx] = norm_names[nid]
alld_idxs_3rd, allc_idxs_3rd, idx_name_3rd

# %%
plt.clf()
fig, ax = plt.subplots(figsize=(6,4))
ax.plot(dat[:,0], dat[:,1], label='Cooperation Level')
for idx,name in idx_name_3rd.items():
  if idx+2 < dat.shape[1]:
    if np.max(dat[:,idx+2]) > 0.05:
      ax.plot(dat[:,0], dat[:,idx+2], label=name)
# sum of ALLD
alld = np.zeros(dat.shape[0])
for idx in alld_idxs_3rd:
  if idx+2 < dat.shape[1]:
    alld += dat[:,idx+2]
allc = np.zeros(dat.shape[0])
for idx in allc_idxs_3rd:
  if idx+2 < dat.shape[1]:
    allc += dat[:,idx+2]
if np.max(alld) > 0.05:
  ax.plot(dat[:,0], alld, label='ALLD')

ax.set_xlim(0,1000)
ax.set_xticks([0,500,1000])
ax.set_xticklabels(['0','500','1000'])
#ax.set_xticklabels(['0','100','200','300','400'])
ax.set_xlabel('time', fontsize=18)
ax.set_ylabel('fractions', fontsize=18)
ax.set_ylim(-0.005,1.005)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
cmap = plt.get_cmap("tab10")
ax.text(990, 0.93, 'Cooperation Level', fontsize=12, color=cmap(0), ha='right', va='top')
ax.text(990, 0.58, r"L1 with $P(B,B)=D$", fontsize=12, color=cmap(1), ha='right', va='bottom')
ax.text(990, 0.40, "L1", fontsize=12, color=cmap(2), ha='right', va='bottom')
ax.text(990, 0.01, "ALLD", fontsize=12, color=cmap(3), ha='right', va='bottom')
#ax.legend()

# %%
fig.savefig('3rd_timeseries.pdf', bbox_inches='tight')
# %%
dat = np.loadtxt('result/2nd_mu0.02_b3_timeseries.dat')

plt.clf()
fig, ax = plt.subplots(figsize=(6,4))
ax.plot(dat[:,0], dat[:,1], label='Cooperation Level')
ax.set_xlabel('time', fontsize=18)
ax.set_ylabel('cooperation level', fontsize=18)
ax.set_ylim(-0.03,1.03)
# %%
norms_2nd = []
with open('2nd_norms.txt', 'r') as file:
  line_number = 0
  for line in file:
    parts = line.strip().split(' ')
    if len(parts) >= 2:
      norms_2nd.append((int(parts[0]), ' '.join(parts[1:])))
norms_2nd

alld_idxs_2nd = []
allc_idxs_2nd = []
idx_name_2nd = {}
for idx,(nid,name) in enumerate(norms_2nd):
  if name == 'AllD':
    alld_idxs_2nd.append(idx)
  if name == 'AllC':
    allc_idxs_2nd.append(idx)
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
    769227: r"L3 with $P(B,B)=C$",
    761035: "L4 BBD",
    859341: r"L7 with $P(B,B)=D$"
  }
  if nid in norm_names:
    idx_name_2nd[idx] = norm_names[nid]
alld_idxs_2nd, allc_idxs_2nd, idx_name_2nd
# %%
plt.clf()
cmap = plt.get_cmap("tab10")
fig, ax = plt.subplots(figsize=(6,4))
ax.plot(dat[:,0], dat[:,1], label='Cooperation Level')
for idx,name in idx_name_2nd.items():
  if idx+2 < dat.shape[1]:
    if np.max(dat[:,idx+2]) > 0.05:
      ax.plot(dat[:,0], dat[:,idx+2], label=name, color=cmap(4))
# sum of ALLD
alld = np.zeros(dat.shape[0])
for idx in alld_idxs_2nd:
  if idx+2 < dat.shape[1]:
    alld += dat[:,idx+2]
allc = np.zeros(dat.shape[0])
for idx in allc_idxs_2nd:
  if idx+2 < dat.shape[1]:
    allc += dat[:,idx+2]

cmap = plt.get_cmap("tab10")
if np.max(alld) > 0.05:
  # using 3rd color in the palette
  ax.plot(dat[:,0], alld, label='ALLD', color=cmap(3))

ax.set_xlim(0,1000)
ax.set_xticks([0,500,1000])
ax.set_xticklabels(['0','500','1000'])
ax.set_xlabel('time', fontsize=18)
ax.set_ylabel('fractions', fontsize=18)
ax.set_ylim(-0.005,1.005)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
cmap = plt.get_cmap("tab10")
ax.text(990, 0.16, 'Cooperation Level', fontsize=12, color=cmap(0), ha='right', va='top')
ax.text(990, 0.82, "ALLD", fontsize=12, color=cmap(3), ha='right', va='bottom')
ax.text(990, 0.195, "L3", fontsize=12, color=cmap(4), ha='right', va='bottom')
#ax.legend()
# %%
fig.savefig('2nd_timeseries.pdf', bbox_inches='tight')
# %%
def plot_2nd_timeseries(ax, datpath):
  dat = np.loadtxt(datpath)
  cmap = plt.get_cmap("tab10")
  ax.plot(dat[:,0], dat[:,1], label='Cooperation Level')
  for idx,name in idx_name_2nd.items():
    if name == 'L3':
      ax.plot(dat[:,0], dat[:,idx+2], label=name, color=cmap(4))
      break
  # sum of ALLD
  alld = np.zeros(dat.shape[0])
  for idx in alld_idxs_2nd:
    if idx+2 < dat.shape[1]:
      alld += dat[:,idx+2]
  ax.plot(dat[:,0], alld, label='ALLD', color=cmap(3))
  ax.set_xlim(0,500)
  #ax.set_xticks([0,100,200,300,400])
  #ax.set_xticklabels(['0','100','200','300','400'])
  ax.set_ylim(-0.005,1.005)

# %%
plt.clf()
fig, axs = plt.subplots(4,4,figsize=(12,12), sharex=True, sharey=True)
for i, mu in enumerate(reversed(['0.01','0.02','0.05','0.1'])):
  for j,b in enumerate(['2','3','4','5']):
    plot_2nd_timeseries(axs[i][j], f'result/2nd_mu{mu}_b{b}_timeseries.dat')
    if i == 0:
      axs[i][j].set_title(r'$b={}$'.format(b), fontsize=18)
    if j == 0:
      axs[i][j].set_ylabel('fractions', fontsize=18)
      axs[i][j].text(-0.4, 0.5, r'$\epsilon={}$'.format(mu), ha='center', va='center', transform=axs[i][j].transAxes, rotation=90, fontsize=18)
    if i == 3:
      axs[i][j].set_xlabel('time', fontsize=18)
    #if i == 0 and j == 3:
      #axs[i][j].legend()
handles, labels = axs[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=12)
plt.tight_layout()
fig.savefig("2nd_all_short_timeseries.pdf", bbox_inches="tight", pad_inches=0.1)
# %%
cmap = plt.get_cmap("tab10")
id_color = {
  765131: 2, #"L1",
  634059: "L2",
  769226: 4, #"L3",
  761034: "L4",
  638154: "L5",
  629962: "L6",
  859333: 6, #"L7", # or 756938
  892101: "L8", # or 625866
  765130: 1, #r"L1 with $P(B,B)=D$",
  765129: "L1 BGD",
  634058: "L2 BBD",
  769227: 5, #r"L3 with $P(B,B)=C$",
  761035: "L4 BBD",
  859341: 7, #r"L7 with $P(B,B)=D$"
}
def plot_3rd_timeseries(ax, datpath):
  dat = np.loadtxt(datpath)
  ax.plot(dat[:,0], dat[:,1], label='Cooperation Level')
  for idx,name in idx_name_3rd.items():
    if idx+2 < dat.shape[1]:
      if np.max(dat[:,idx+2]) > 0.05:
        ax.plot(dat[:,0], dat[:,idx+2], label=name, color=cmap(id_color[norms_3rd[idx][0]]) )
  # sum of ALLD
  alld = np.zeros(dat.shape[0])
  for idx in alld_idxs_3rd:
    if idx+2 < dat.shape[1]:
      alld += dat[:,idx+2]
  if np.max(alld) > 0.05:
    ax.plot(dat[:,0], alld, label='ALLD', color=cmap(3))
  ax.set_ylim(-0.005,1.005)
  #ax.legend()

# %%
fig, axs = plt.subplots(4,4,figsize=(12,12), sharex=False, sharey=True)
for i, mu in enumerate(reversed(['0.01','0.02','0.05','0.1'])):
  for j,b in enumerate(['2','3','4','5']):
    plot_3rd_timeseries(axs[i][j], f'result/3rd_mu{mu}_b{b}_timeseries.dat')
    if i == 0:
      axs[i][j].set_title(r'$b={}$'.format(b), fontsize=18)
    if i == 2 or i == 3:
      axs[i][j].set_xlim(0,450)
      axs[i][j].set_xticks([0,100,200,300,400])
      axs[i][j].set_xticklabels(['0','100','200','300','400'])
    elif i == 0 and j == 0:
      axs[i][j].set_xlim(0,18000)
      axs[i][j].set_xticks([0,5000,10000,15000])
      axs[i][j].set_xticklabels(['0','5000','10000','15000'])
    else:
      axs[i][j].set_xlim(0,1800)
      axs[i][j].set_xticks([0,500,1000,1500])
      axs[i][j].set_xticklabels(['0','500','1000','1500'])
    if j == 0:
      axs[i][j].set_ylabel('fractions', fontsize=18)
      axs[i][j].text(-0.4, 0.5, r'$\epsilon={}$'.format(mu), ha='center', va='center', transform=axs[i][j].transAxes, rotation=90, fontsize=18)
    if i == 3:
      axs[i][j].set_xlabel('time', fontsize=18)
    #if i == 0 and j == 3:
      #axs[i][j].legend()
#handles, labels = axs[0, 0].get_legend_handles_labels()
#fig.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=12)
handles, labels = [], []
for ax in fig.axes:
  for handle, label in zip(*ax.get_legend_handles_labels()):
    if label not in labels:  # To avoid duplicate labels
      handles.append(handle)
      labels.append(label)
# rearrange the order
handles[1], handles[2], handles[3] = handles[2], handles[3], handles[1]
labels[1], labels[2], labels[3] = labels[2], labels[3], labels[1]
fig.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=12)
plt.tight_layout()
fig.savefig("3rd_all_short_timeseries.pdf", bbox_inches="tight", pad_inches=0.2)


# %%
