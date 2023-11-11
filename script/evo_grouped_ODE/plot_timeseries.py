#%%
import matplotlib.pyplot as plt
import numpy as np
# %%
# Load data
dat = np.loadtxt('timeseries.dat')

plt.clf()
fig, ax = plt.subplots(figsize=(6,4))
ax.plot(dat[:,0], dat[:,1], label='Cooperation Level')
ax.set_xlabel('time', fontsize=18)
ax.set_ylabel('cooperation level', fontsize=18)
ax.set_ylim(-0.03,1.03)

# %%
fig.savefig('timeseries.png')
# %%
# import msgpack
# 
# norm_ids = msgpack.load(open('../fix_prob_results/third_order_mu0.01/fixation_probs_15.msgpack', 'rb'))['norm_ids']
# norm_names = {
#   765131: "L1",
#   634059: "L2",
#   769226: "L3",
#   761034: "L4",
#   638154: "L5",
#   629962: "L6",
#   859333: "L7", # or 756938
#   892101: "L8", # or 625866
#   765130: "L1 BBD",
#   765129: "L1 BGD",
#   634058: "L2 BBD",
#   769227: "L3 BBC"
# }
# 
# name_index = {}
# for nid,name in norm_names.items():
#   idx = norm_ids.index(nid)
#   name_index[name] = idx
# name_index

# %%
name_index = {
  'L1': 1107,
  'L2': 751,
  'L3': 1122,
  'L4': 1090,
  'L5': 766,
  #'L6': 735,
  'L7': 1383,
  #'L8': 1489,
  'L1 BBD': 1106,
  'L1 BGD': 1105,
  'L2 BBD': 750,
  'L3 BBC': 1123
}

# %%
plt.clf()
fig, ax = plt.subplots(figsize=(6,4))
ax.plot(dat[:,0], dat[:,1], label='Cooperation Level')
for name,idx in name_index.items():
  ax.plot(dat[:,0], dat[:,idx+2], label=name)
ax.set_xlabel('time', fontsize=18)
ax.set_ylabel('cooperation level', fontsize=18)
ax.set_ylim(-0.03,1.03)
ax.legend()

# %%
fig.savefig('timeseries_with_norms.png')
# %%
