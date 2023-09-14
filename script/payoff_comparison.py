# %%
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import os

# %%
script_path = os.path.dirname(os.path.abspath(__file__))
exe_path = os.path.join(script_path, '..', 'cmake-build-release', 'inspect_EvolPrivRepGame')
print(exe_path)

# %%
benefit = 5

def calc_payoffs(resident, mutant):
  #t_measure = 1000
  t_measure = 10000
  out = subprocess.check_output([exe_path, '-j', f'{{"N":50,"q":1,"t_measure":{t_measure},"benefit":{benefit}}}', resident, mutant], universal_newlines=True)

  dat = np.loadtxt(out.splitlines(), delimiter=' ', skiprows=0)
  return dat

# %%
norms = ['L1','L2','L3','L4','L5','L6','L7','L8','AllC','AllD']
#norms = ['L1','AllC','AllD']
num_norms = len(norms)
dat_array = [[None] * num_norms for i in range(num_norms)]
# take all combinations of norms
for i in range(num_norms):
  print(i)
  for j in range(i+1, num_norms):
    dat = calc_payoffs(norms[i], norms[j])
    dat_array[i][j] = dat

# %%
dat_array
# %%
# make plots
fig, axs = plt.subplots(num_norms, num_norms, figsize=(4*num_norms,4*num_norms))

def make_subplots(dat, ax, ax2, resident, mutant):
  ax.plot(dat[:, 0], dat[:, 1], label=f"resident: {resident}")
  ax.plot(dat[:, 0], dat[:, 2], label=f"mutant: {mutant}")
  ax.set_xlabel('# of mutants')
  ax.set_ylabel('payoff')
  ax.set_ylim(-1,benefit)
  ax.legend()

  ax2.plot(dat[:, 0], dat[::-1, 2], label=f"resident: {mutant}")
  ax2.plot(dat[:, 0], dat[::-1, 1], label=f"mutant: {resident}")
  ax2.set_xlabel('# of mutants')
  ax2.set_ylabel('payoff')
  ax2.set_ylim(-1,5)
  ax2.legend()

for i in range(num_norms):
  print(i)
  fig.delaxes(axs[i, i])
  for j in range(i+1, num_norms):
    ax = axs[i, j]
    ax2 = axs[j, i]
    make_subplots(dat_array[i][j], ax, ax2, norms[i], norms[j])

# %%
fig.savefig(f'payoff_comparison_b{benefit}.pdf')
# %%
