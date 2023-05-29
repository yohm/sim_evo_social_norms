#%%
import matplotlib.pyplot as plt
import numpy as np
# %%
# Load data
data = np.loadtxt('timeseries.dat', delimiter=' ')
data
# %%
# Plot data
plt.clf()
plt.plot(data[:,0], data[:,1], label='Cooperation Level')
plt.plot(data[:,0], data[:,2], label='L1')
plt.plot(data[:,0], data[:,3], label='L1v')
plt.plot(data[:,0], data[:,4], label='L3')
plt.xlabel('Time')
plt.ylabel('Cooperation Level')
plt.ylim(0,1)
plt.legend()
plt.savefig('timeseries.png')
# %%
