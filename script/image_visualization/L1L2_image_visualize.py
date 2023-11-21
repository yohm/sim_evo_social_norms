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
  out = subprocess.check_output([exe_path, '-j', f'{{"N":50,"q":1,"mu_assess1":0.01,"t_measure":{t_measure},"benefit":{benefit}}}', resident, mutant], universal_newlines=True)

  dat = np.loadtxt(out.splitlines(), delimiter=' ', skiprows=0)
  return dat

# %%
dat = calc_payoffs('L1', 'L2')
# %%
# Read the text file
with open('L1L2_image.txt', 'r') as file:
    lines = file.readlines()

# %%
# Map characters '.' and 'x' to colors
color_map = {
    '.': [1.0, 1.0, 1.0],  # White color for '.'
    'x': [0.2, 0.2, 0.2]   # Black color for 'x'
}

# Convert characters to colors and create a 2D numpy array
pixels = np.array([[color_map[char] for char in line.strip()] for line in lines])
pixels

# %%
# Display the image using matplotlib
plt.clf()
fig, ax = plt.subplots(figsize=(6, 4))

ax.plot(dat[:, 0], dat[:, 1], label=f"L1")
ax.plot(dat[:, 0], dat[:, 2], label=f"L2")
ax.set_xlabel('# of L2 players', fontsize=18)
ax.set_ylabel('payoffs', fontsize=18)
ax.set_ylim(-1,benefit)
ax.legend()
#ax.text(-0.12, 1.12, 'a', transform=ax.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

fig.savefig('L1L2_payoff.pdf', bbox_inches='tight', pad_inches=0.2)

# %%
plt.clf()
fig, ax1 = plt.subplots(figsize=(4, 4))

ax1.imshow(pixels)
ax1.axis('on')  # Turn off axis for better visualization

# 1Draw a rectangle to indicate the region
ax1.axvline(x=24.5, color='y', linewidth=1, linestyle='--')
ax1.axhline(y=24.5, color='y', linewidth=1, linestyle='--')
ax1.text(12, -3, 'L1', fontsize=16, color='black', ha='center', va='center')
ax1.text(37, -3, 'L2', fontsize=16, color='black', ha='center', va='center')
ax1.text(-3, 12, 'L1', fontsize=16, color='black', ha='center', va='center', rotation=90)
ax1.text(-3, 37, 'L2', fontsize=16, color='black', ha='center', va='center', rotation=90)

ax1.text(20, 2.4, "95.4%", fontsize=12, bbox=dict(facecolor='white', alpha=0.8, edgecolor='black', boxstyle='round,pad=0.3'), ha='center', va='center')
ax1.text(45, 2.4, "85.9%", fontsize=12, bbox=dict(facecolor='white', alpha=0.9, edgecolor='black', boxstyle='round,pad=0.3'), ha='center', va='center')
ax1.text(20, 27.4, "84.6%", fontsize=12, bbox=dict(facecolor='white', alpha=0.9, edgecolor='black', boxstyle='round,pad=0.3'), ha='center', va='center')
ax1.text(45, 27.4, "80.0%", fontsize=12, bbox=dict(facecolor='white', alpha=0.9, edgecolor='black', boxstyle='round,pad=0.3'), ha='center', va='center')
#ax1.text(-0.12, 1.12, 'b', transform=ax1.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

# ticks off
ax1.set_xticks([])
ax1.set_yticks([])
# %%
fig.savefig('L1L2_image.pdf', bbox_inches='tight', pad_inches=0.2)

# %%
