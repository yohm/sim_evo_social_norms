# %%
import numpy as np
import matplotlib.pyplot as plt
import os,subprocess
import json

# %%
script_path = os.path.dirname(os.path.abspath(__file__))
exe_path = os.path.join(script_path, '..', 'cmake-build-release', 'inspect_EvolPrivRepGame')

# %%
def calc_payoffs(mu_assess=0.01, mu_impl=0.0, t_measure=10000, t_init=10000, q=1.0, n=50, benefit=5.0, sigma_in=1.0, resident='L1'):
  out = subprocess.check_output([exe_path, '-j', f'{{"N":{n},"mu_assess1":{mu_assess},"mu_assess2":0.0,"mu_impl":{mu_impl},"mu_percept":0.0,"q":{q},"t_init":{t_init},"t_measure":{t_measure},"benefit":{benefit},"beta":{sigma_in}}}', resident], universal_newlines=True)
  return json.loads(out)

# %%
j = calc_payoffs(n=5, benefit=5.0, sigma_in=1.0, resident='L1')
j
# %%

mu_list = [0.01, 0.03, 0.05, 0.1]
benefit_list = [1.5,2.0,2.5,3.0,4.0,5.0]

N = 50
resident = 'L1'
dats = {}
for mu_i,mu in enumerate(mu_list):
  dats[mu_i] = []
  for b_i,benefit in enumerate(benefit_list):
    out = calc_payoffs(mu_assess=mu, benefit=benefit, n=N, resident=resident)
    dats[mu_i].append(out['eq_cooperation_level'])

# %%
dats

# %%
plt.clf()
color_map = plt.get_cmap('viridis')
plt.xlabel('benefit')
plt.ylabel('cooperation level')
plt.xlim([1.2,5.2])
plt.ylim([0.0,1.0])
for mu_i,mu in enumerate(mu_list):
  plt.plot(benefit_list, dats[mu_i], label=f'$\mu_a={mu}$', marker='o', color=color_map(mu_i/len(mu_list)))
plt.legend(loc='upper left')
plt.savefig(f'three_species_{resident}.pdf')


# %%
# %%
mu_list

# %%
