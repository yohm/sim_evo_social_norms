# %%
import msgpack
# %%

# load resluts file in msgpack format
# find msgpack files in the current directory
import glob
msgpack_files = glob.glob('**/*.msgpack')
print(msgpack_files)

# %%
map_to_param = {}
for msgpack_file in msgpack_files:
    with open(msgpack_file, 'rb') as f:
        data = msgpack.unpack(f)
    print(data.keys())
    params = data['evoparam']
    params['benefit'] = data['benefit']
    params['sigma_in'] = data['sigma_in']
    params['norm_set'] = data['norm_set']
    map_to_param[msgpack_file] = params
map_to_param
# %%
# print in json
# write to file 'path_to_param.json'
import json
open('path_to_param.json', 'w').write(json.dumps(map_to_param, indent=2))


# %%
