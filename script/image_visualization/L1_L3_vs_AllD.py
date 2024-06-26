# %%
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import subprocess
import os

# %%
script_path = os.path.dirname(os.path.abspath(__file__))

# %%
benefit = 5
t_measure = 100000

def calc_payoffs(resident, mutant):
  exe_path = os.path.join(script_path, '..', '..', 'cmake-build-release', 'inspect_EvolPrivRepGame')
  out = subprocess.check_output([exe_path, '-j', f'{{"N":50,"q":1,"mu_assess1":0.02,"t_measure":{t_measure},"benefit":{benefit}}}', resident, mutant], universal_newlines=True)

  dat = np.loadtxt(out.splitlines(), delimiter=' ', skiprows=0)
  return dat

# %%
def calc_image(resident, mutant):
  exe_path = os.path.join(script_path, '..', '..', 'cmake-build-release', 'inspect_PrivRepGame')
  subprocess.run([exe_path, '-j', f'{{"q":1,"mu_assess1":0.02,"t_measure":10000}}', resident, "25", mutant, "25", '-g'], universal_newlines=True)
  with open('image.txt', 'r') as file:
    lines = file.readlines()
  color_map = {
    '.': [1.0, 1.0, 1.0],  # White color for '.'
    'x': [0.2, 0.2, 0.2]   # Black color for 'x'
  }
  pixels = np.array([[color_map[char] for char in line.strip()] for line in lines])
  return pixels


# %%
dat1 = calc_payoffs('L1', 'AllD')
dat2 = calc_payoffs('L3', 'AllD')

# %%
pixel1 = calc_image('L1', 'AllD')
pixel2 = calc_image('L3', 'AllD')
pixel1,pixel2

# %%
plt.clf()
fig, ((ax2,ax3),(ax0,ax1)) = plt.subplots(figsize=(7, 6.3), nrows=2, ncols=2)

ax0.imshow(pixel1[0:25])
ax0.set_xticks([])
ax0.set_yticks([])
ax0.axvline(x=24.5, color='y', linewidth=2, linestyle='--')
ax0.text(-3, 12, "L1's opinion", fontsize=12, color='black', ha='center', va='center', rotation=90)
ax0.text(12, -3, 'L1', fontsize=12, color='black', ha='center', va='center')
ax0.text(37, -3, 'ALLD', fontsize=12, color='black', ha='center', va='center')
ax0.set_ylim(24.5,-0.5)
#ax0.text(-0.08, 1.3, 'c', transform=ax0.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

ax1.imshow(pixel2[0:25])
ax1.set_xticks([])
ax1.set_yticks([])
ax1.axvline(x=24.5, color='y', linewidth=2, linestyle='--')
ax1.text(-3, 12, "L3's opinion", fontsize=12, color='black', ha='center', va='center', rotation=90)
ax1.text(12, -3, 'L3', fontsize=12, color='black', ha='center', va='center')
ax1.text(37, -3, 'ALLD', fontsize=12, color='black', ha='center', va='center')
ax1.set_ylim(24.5,-0.5)
#ax1.text(-0.08, 1.3, 'd', transform=ax1.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

cmap = plt.get_cmap('tab10')
ax2.plot(dat1[:, 0], dat1[:, 1], label=f"L1")
ax2.plot(dat1[:, 0], dat1[:, 2], label=f"ALLD")
ax2.set_xlabel('# of ALLD players', fontsize=12)
ax2.set_ylabel('payoffs', fontsize=14)
ax2.set_ylim(-1,benefit)
ax2.set_xlim(0,50)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.text(15, 2.75, 'L1', fontsize=12, color=cmap(0), ha='center', va='bottom')
ax2.text(15, 0.2, 'ALLD', fontsize=12, color=cmap(1), ha='center', va='bottom')
ax2.tick_params(axis='both', which='major', labelsize=14)
ax2.set_xticks([0, 25, 50])
#ax2.legend()
#ax2.text(-0.08, 1.16, 'a', transform=ax2.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')


ax3.plot(dat2[:, 0], dat2[:, 1], label=f"L3")
ax3.plot(dat2[:, 0], dat2[:, 2], label=f"ALLD")
ax3.set_xlabel('# of ALLD players', fontsize=12)
ax3.set_ylim(-1,benefit)
ax3.set_xlim(0,50)
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.set_yticklabels([])
ax3.text(15, 2.75, 'L1', fontsize=12, color=cmap(0), ha='center', va='bottom')
ax3.text(15, 1.0, 'ALLD', fontsize=12, color=cmap(1), ha='center', va='bottom')
ax3.tick_params(axis='both', which='major', labelsize=14)
ax3.set_xticks([0, 25, 50])
#ax3.legend()
#ax3.text(-0.08, 1.16, 'b', transform=ax3.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

# %%
fig.savefig('L1_L3_vs_AllD.pdf', bbox_inches='tight', pad_inches=0.2)
# %%
plt.clf()
fig, ax = plt.subplots(figsize=(6, 4))

ax.imshow(pixel1[0:25])
ax.set_xticks([])
ax.set_yticks([])
ax.axvline(x=24.5, color='y', linewidth=2, linestyle='--')
ax.text(-3, 12, "L1's opinion", fontsize=18, color='black', ha='center', va='center', rotation=90)
ax.text(12, -3, 'L1', fontsize=18, color='black', ha='center', va='center')
ax.text(37, -3, 'ALLD', fontsize=18, color='black', ha='center', va='center')
ax.set_ylim(24.5,-0.5)

# %%
fig.savefig('L1_vs_AllD_image.pdf', bbox_inches='tight')

# %%
plt.clf()
fig, ax = plt.subplots(figsize=(6, 4))

ax.imshow(pixel2[0:25])
ax.set_xticks([])
ax.set_yticks([])
ax.axvline(x=24.5, color='y', linewidth=2, linestyle='--')
ax.text(-3, 12, "L3's opinion", fontsize=18, color='black', ha='center', va='center', rotation=90)
ax.text(12, -3, 'L3', fontsize=18, color='black', ha='center', va='center')
ax.text(37, -3, 'ALLD', fontsize=18, color='black', ha='center', va='center')
ax.set_ylim(24.5,-0.5)
# %%
fig.savefig('L3_vs_AllD_image.pdf', bbox_inches='tight')
# %%
plt.clf()
fig, ax = plt.subplots(figsize=(6, 4))

cmap = plt.get_cmap('tab10')
ax.plot(dat1[:, 0], dat1[:, 1], label=f"L1")
ax.plot(dat1[:, 0], dat1[:, 2], label=f"ALLD")
ax.set_xlabel('# of ALLD players', fontsize=18)
ax.set_ylabel('payoffs', fontsize=18)
ax.set_ylim(-1,benefit)
ax.set_xlim(0,50)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.text(15, 2.75, 'L1', fontsize=18, color=cmap(0), ha='center', va='bottom')
ax.text(15, 0.2, 'ALLD', fontsize=18, color=cmap(1), ha='center', va='bottom')
ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_xticks([0, 25, 50])
#ax2.legend()
# %%
fig.savefig('L1_vs_AllD_payoff.pdf', bbox_inches='tight')

# %%
plt.clf()
fig, ax = plt.subplots(figsize=(6, 4))

ax.plot(dat2[:, 0], dat2[:, 1], label=f"L3")
ax.plot(dat2[:, 0], dat2[:, 2], label=f"ALLD")
ax.set_xlabel('# of ALLD players', fontsize=18)
ax.set_ylabel('payoffs', fontsize=18)
ax.set_ylim(-1,benefit)
ax.set_xlim(0,50)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.text(15, 2.75, 'L3', fontsize=18, color=cmap(0), ha='center', va='bottom')
ax.text(15, 1.0, 'ALLD', fontsize=18, color=cmap(1), ha='center', va='bottom')
ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_xticks([0, 25, 50])
# %%
fig.savefig("L3_vs_AllD_payoff.pdf", bbox_inches="tight")
# %%
