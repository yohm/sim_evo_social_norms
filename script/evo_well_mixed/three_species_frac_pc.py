# %%
import numpy as np
import matplotlib.pyplot as plt

# %%
# run in advance
# ../../cmake-build-release/main_well_mixed_evo -3 ../fix_prob_results/third_order_mu0.01/fixation_probs_15.msgpack > frac_pc.dat
dat = np.loadtxt('frac_pc.dat')

# %%
leading_eight = [
    765131, # L1
    634059, # L2
    769226, # L3
    761034, # L4
    638154, # L5
    629962, # L6
    859333, # L7 # or 756938
    892101, # L8 # or 625866
]
secondary_sixteen = [
    957633, # S1
    695489, # S2
    617673, # S3
    621768, # S4
    826561, # S5
    646344, # S6
    650441, # S7
    654536, # S8
    924865, # S9
    744648, # S10
    748745, # S11
    752840, # S12
    793793, # S13
    777416, # S14
    781513, # S15
    785608, # S16
]

# exclude data points whose dat[0] is in the leading eight
non_dat = dat[np.isin(dat[:,0], leading_eight+secondary_sixteen, invert=True)]
l8_dat = dat[np.isin(dat[:,0], leading_eight)]
s16_dat = dat[np.isin(dat[:,0], secondary_sixteen)]
non_dat.shape, l8_dat.shape, s16_dat.shape


# %%
plt.clf()
fig, ax = plt.subplots(figsize=(6, 4), ncols=1)
ax.plot(non_dat[:,2], non_dat[:,1], '.', color='#6883BA')
ax.set_xlim(-0.02,1.02)
ax.set_ylim(-0.02,1.02)
ax.set_xlabel('self-cooperation level', fontsize=12)
ax.set_ylabel('equilibrium fraction', fontsize=12)

ax.plot(s16_dat[:,2], s16_dat[:,1], 's', color='#3D3B8E')
ax.plot(l8_dat[:,2], l8_dat[:,1], 'o', color='#E072A4')

ax.text(0.97, 0.76, 'L1', fontsize=12, color='#222222', horizontalalignment='right', verticalalignment='bottom')
ax.text(0.86, 0.97, 'L2', fontsize=12, color='#222222', horizontalalignment='left', verticalalignment='top')
ax.text(0.97, 0.07, 'L3', fontsize=12, color='#222222', horizontalalignment='right', verticalalignment='bottom')
ax.text(0.97, 0.05, 'L4', fontsize=12, color='#222222', horizontalalignment='right', verticalalignment='center')
ax.text(0.87, 0.15, 'L5', fontsize=12, color='#222222', horizontalalignment='right', verticalalignment='bottom')
ax.text(0.49, 0.01, 'L6', fontsize=12, color='#222222', horizontalalignment='right', verticalalignment='bottom')
ax.text(0.97, 0.61, 'L7', fontsize=12, color='#222222', horizontalalignment='right', verticalalignment='bottom')
ax.text(0.10, 0.09, 'L8', fontsize=12, color='#222222', horizontalalignment='right', verticalalignment='bottom')

# %%
# sort l8_dat according to the order of the leading eight
l8_dat
# %%
fig.savefig('three_species_frac_pc.pdf', bbox_inches='tight')

# %%
