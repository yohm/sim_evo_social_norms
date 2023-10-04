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

fig.savefig('timeseries.png')
# %%
