# %%
import numpy as np
import matplotlib.pyplot as plt

# %%
plt.clf()
fig,ax = plt.subplots(1,1,figsize=(6,4))

ax.set_xlabel("# of groups", fontsize=18)
ax.set_ylabel("fraction", fontsize=18)

ax.set_xlim([-1, 102])
ax.set_ylim([-0.02, 1.02])
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

x = [1, 2, 5, 10, 20, 50, 100]
y = [0.13629008688024002, 0.23366875558677985, 0.7575360129302602, 0.857491606689728, 0.9006166220448348, 0.9327120089287891, 0.9434438805764293]
yerr = [0.004130113610016797, 0.003417787458878201, 0.01130636348408964, 0.008558765204420475, 0.0027257720734248213, 0.0031557098145661502, 0.0020613329477889857]

l1_frac = [0.00833864, 0.02928902, 0.237683352, 0.32978575600000004, 0.5102484985, 0.6971432824, 0.8756457400000001]
l1_frac_err = [0.0013999121702767245, 0.0032062611493780587, 0.015024817499202158, 0.02475288612548424, 0.025764470439219856, 0.03760777675117192, 0.019020730337499396]

ax.errorbar(x, y, yerr=yerr, marker='o', label="cooperation level")
ax.errorbar(x, l1_frac, yerr=l1_frac_err, marker='o', label="L1 fraction")

cmap = plt.get_cmap("tab10")
ax.text(95, 0.94, 'Cooperation Level', fontsize=12, color=cmap(0), ha='right', va='bottom')
ax.text(95, 0.77, 'L1 fraction', fontsize=12, color=cmap(1), ha='right', va='top')

# %%
fig.savefig("pc_M_dep.pdf", bbox_inches="tight")

# %%
plt.clf()
fig,ax = plt.subplots(1,1,figsize=(6,4))

ax.set_xlabel("time", fontsize=18)
ax.set_ylabel("fraction", fontsize=18)

ax.set_ylim([-0.02, 1.02])
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# load file from timeseries_M100.dat
dat = np.loadtxt("timeseries_M100.dat")

ax.plot(dat[:,0], dat[:,1], label="cooperation level")
ax.plot(dat[:,0], dat[:,2], label="L1 fraction")

ax.legend(loc='lower left', fontsize=12)
# %%
fig.savefig("timeseries_M100.pdf", bbox_inches="tight")
# %%
