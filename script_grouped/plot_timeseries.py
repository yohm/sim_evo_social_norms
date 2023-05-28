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
plt.plot(data[:,0], data[:,1])
plt.xlabel('Time')
plt.ylabel('Cooperation Level')
plt.savefig('timeseries.pdf')
# %%
