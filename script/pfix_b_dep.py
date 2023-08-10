# %%
import numpy as np
import matplotlib.pyplot as plt
import os,subprocess

# %%
script_path = os.path.dirname(os.path.abspath(__file__))
exe_path = os.path.join(script_path, '..', 'cmake-build-release', 'main_fix_probs_param_dep')

# %%
def calc_payoffs(mu_assess=0.01, n_steps=10000, q=1.0, n=50, resident='AllD', mutant='L1'):
  out = subprocess.check_output([exe_path, '-j', f'{{"N":{n},"q":{q},"n_steps":{n_steps},"mu_assess":{mu_assess},"n_steps":{n_steps}}}', resident, mutant], universal_newlines=True)
  dat = np.loadtxt(out.splitlines())
  return dat

# %%
# %%
def make_subplot(ax, dat, ylim_max=0.1):
  y = dat[:, 1]
  unique_y = np.unique(dat[:, 1])
  cmap = plt.cm.get_cmap('viridis', len(unique_y)).reversed()
  for i,y_val in enumerate(unique_y):
    # Filter data for the current y_val
    filtered_data = dat[dat[:, 1] == y_val]

    x = filtered_data[:, 0]  # X values (first column)
    z = filtered_data[:, 2]  # Y values (third column)
    ax.plot(x, z, linestyle='-', label=f'$\sigma = {y_val}$', color=cmap(i))
  ax.set_xlim(1, 6)
  ax.set_ylim(0, ylim_max)
  #ax.legend()
  ax.grid(True)
# %%

mu_list = [0.01, 0.03, 0.05, 0.1]
q_list = [1.0, 0.9, 0.5, 0.1]

N = 30
resident = 'AllD'
mutant = 'L1'
dats = {}
for mu_i,mu in enumerate(mu_list):
  for q_i,q in enumerate(q_list):
    dat = calc_payoffs(mu_assess=mu, q=q, n=N, resident=resident, mutant=mutant)
    dats[(mu_i, q_i)] = dat

# %%
fig,axs = plt.subplots(len(mu_list), len(q_list), figsize=(12,12), sharex=True, sharey=True)

for mu_i,mu in enumerate(mu_list):
  for q_i,q in enumerate(q_list):
    dat = dats[(mu_i, q_i)]
    ax = axs[mu_i,q_i]
    #dat = calc_payoffs(mu_assess=mu, q=q)
    ylim_max = 5 * 1.0 / N
    make_subplot(ax, dat, ylim_max=ylim_max)
    if mu_i == 0:
      ax.set_title(f'$q = {q}$')
    if mu_i == len(mu_list)-1:
      ax.set_xlabel('benefit')
    if q_i == 0:
      title_text = f"$\mu = {mu}$"
      ax.text(-0.4, 0.5, title_text, transform=ax.transAxes,
        rotation=90, va='center', ha='center',
        fontsize=12, fontweight='bold')
      ax.set_ylabel('fixation probability')

# %%
fig.savefig(f"pfix_b_dep_N{N}_{resident}.pdf")