# %%
import json
import os

# %%

# load file from fix_prob_results/path_to_param.json
json_path = os.path.dirname(__file__) + '/../fix_prob_results/path_to_param.json'
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
comp_keys = ['norm_set','N','benefit','sigma_in','mu_assess1','mu_impl','q']
p1 = [input_json[k] for k in comp_keys]
fix_prob_path = ""
for k,v in fix_param.items():
    p2 = [v[k] for k in comp_keys]
    if p1 == p2:
        fix_prob_path = f"{os.path.dirname(__file__)}/fix_prob_results/{k}"
        break
p1,p2,fix_prob_path
# %%
if fix_prob_path == "":
    print("No matching parameter found")
    exit(1)
# %%
sigma_out, mut_r, T_max, dt = input_json['sigma_out'], input_json['mut_r'], input_json['T_max'], input_json['dt']
print(f"{fix_prob_path} {sigma_out} {mut_r} {T_max} {dt}")
# %%
