# %%
import json
import os

# %%

# load file from fix_prob_results/path_to_param.json
fix_prob_results_dir = os.path.dirname(__file__) + '/../fix_prob_results/'
json_path = fix_prob_results_dir + 'path_to_param.json'
json_path
#%%
with open(json_path) as f:
    fix_param = json.load(f)
fix_param
# %%
with open('_input.json') as f:
    input_json = json.load(f)
input_json
# %%
comp_keys = ['norm_set','N','benefit','sigma_in','mu_percept', 'mu_assess1','q']
p1 = [input_json[k] for k in comp_keys]
fix_prob_path = ""
for k,v in fix_param.items():
    p2 = [v[k] for k in comp_keys]
    if p1 == p2:
        fix_prob_path = f"{fix_prob_results_dir}/{k}"
        break
p1,p2,fix_prob_path
# %%
if fix_prob_path == "":
    print("No matching parameter found")
    exit(1)
# %%
arg_json_keys = ['benefit','M','T_init','T_measure','mut_r','sigma_out','_seed']
d = {k:input_json[k] for k in arg_json_keys}
d['seed'] = d.pop('_seed')  # change key name
d

# %%
print(f"{fix_prob_path} -j '{json.dumps(d)}'")
# %%
