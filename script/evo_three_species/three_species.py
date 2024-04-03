# %%
import numpy as np
import matplotlib.pyplot as plt
import os,subprocess
import json,pickle

# %%
script_path = os.path.dirname(os.path.abspath(__file__))
exe_path = os.path.join(script_path, '..', '..', 'cmake-build-release', 'inspect_EvolPrivRepGame')

# %%
def calc_payoffs(mu_assess=0.01, mu_impl=0.0, t_measure=10000, t_init=10000, q=1.0, n=50, benefit=5.0, sigma_in=1.0, resident='L1'):
  out = subprocess.check_output([exe_path, '-j', f'{{"N":{n},"mu_assess1":{mu_assess},"mu_assess2":0.0,"mu_impl":{mu_impl},"mu_percept":0.0,"q":{q},"t_init":{t_init},"t_measure":{t_measure},"benefit":{benefit},"beta":{sigma_in}}}', resident], universal_newlines=True)
  return json.loads(out)

# %%
# j = calc_payoffs(n=5, benefit=5.0, sigma_in=1.0, resident='L2')
# j

# %%
# if thre is a file 'three_species.pickle', load it
# otherwise, calculate and save it
dat = {}
if os.path.exists('three_species.pickle'):
  with open('three_species.pickle', 'rb') as f:
    dat = pickle.load(f)
dat

# %%
mu_list = [0.01, 0.02, 0.05, 0.1]
benefit_list = [1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0]

N = 50
leading_eight = ['L1','L2','L3','L4','L5','L6','L7','L8']
secondary_sixteen = ['S1','S2','S3','S4','S5','S6','S7','S8','S9','S10','S11','S12','S13','S14','S15','S16']
updated = False
for resident in (leading_eight + secondary_sixteen):
  for mu_i,mu in enumerate(mu_list):
    key = (N,mu,resident)
    if key in dat:
      continue
    dat[key] = []
    updated = True
    for b_i,benefit in enumerate(benefit_list):
      out = calc_payoffs(mu_assess=mu, benefit=benefit, n=N, resident=resident)
      dat[key].append(out['eq_cooperation_level'])
dat

# %%
# save dat to a file
if updated:
  with open('three_species.pickle', 'wb') as f:
    pickle.dump(dat, f)


# %%
plt.clf()
resident = 'L1'
color_map = plt.get_cmap('viridis')
plt.xlabel('benefit')
plt.ylabel('cooperation level')
plt.xlim([1.2,5.2])
plt.ylim([0.0,1.0])
for mu_i,mu in enumerate(mu_list):
  key = (N,mu,resident)
  plt.plot(benefit_list, dat[key], label=f'$\mu_a={mu}$', marker='o', color=color_map(mu_i/len(mu_list)))
plt.legend(loc='upper left')


# %%
fig,axs = plt.subplots(2, 4, figsize=(24,12), sharex=True, sharey=True)
plt.subplots_adjust(wspace=0.02, hspace=0.2)
custom_x_ticks = [2.0, 3.0, 4.0, 5.0]

for i, (ax, resident) in enumerate(zip(axs.flat, leading_eight)):
  ax.set_title(resident, fontsize=26)
  for mu_i,mu in enumerate(mu_list):
    key = (N,mu,resident)
    ax.plot(benefit_list, dat[key], label=f'$\epsilon={mu}$', marker='o', color=color_map(mu_i/len(mu_list)))
    ax.set_xlim([1.2,5.2])
    ax.set_xticks(custom_x_ticks)
    ax.set_xticklabels(custom_x_ticks, fontsize=16)
    ax.set_yticks([0.0,0.2,0.4,0.6,0.8,1.0])
    ax.set_yticklabels([0.0,0.2,0.4,0.6,0.8,1.0], fontsize=16)
    ax.set_ylim([-0.02,1.02])
    if i == 3:
      ax.legend(loc='upper right', fontsize=20)
    if i >= 4:
      ax.set_xlabel('benefit', fontsize=24)
    if i % 4 == 0:
      ax.set_ylabel('cooperation level', fontsize=24)

# %%
fig.savefig('three_species_L8.pdf', bbox_inches='tight', pad_inches=0.3)
# %%

fig,axs = plt.subplots(4, 4, figsize=(24,24), sharex=True, sharey=True)
plt.subplots_adjust(wspace=0.02, hspace=0.2)
custom_x_ticks = [2.0, 3.0, 4.0, 5.0]

for i, (ax, resident) in enumerate(zip(axs.flat, secondary_sixteen)):
  ax.set_title(resident, fontsize=26)
  for mu_i,mu in enumerate(mu_list):
    key = (N,mu,resident)
    ax.plot(benefit_list, dat[key], label=f'$\epsilon={mu}$', marker='o', color=color_map(mu_i/len(mu_list)))
    ax.set_xlim([1.2,5.2])
    ax.set_xticks(custom_x_ticks)
    ax.set_xticklabels(custom_x_ticks, fontsize=16)
    ax.set_yticks([0.0,0.2,0.4,0.6,0.8,1.0])
    ax.set_yticklabels([0.0,0.2,0.4,0.6,0.8,1.0], fontsize=16)
    ax.set_ylim([-0.02,1.02])
    if i == 3:
      ax.legend(loc='upper right', fontsize=20)
    if i >= 12:
      ax.set_xlabel('benefit', fontsize=24)
    if i % 4 == 0:
      ax.set_ylabel('cooperation level', fontsize=24)
# %%
fig.savefig('three_species_S16.pdf', bbox_inches='tight', pad_inches=0.3)

# %%
