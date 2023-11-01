# %%
import numpy as np
import matplotlib.pyplot as plt
import os,sys
from matplotlib.cm import ScalarMappable

# %%
script_path = os.path.dirname(os.path.abspath(__file__))
if 'ipykernel' in sys.modules:
    # we are in a jupyter notebook
    input_file = os.path.join(script_path, '../../script_grouped/fix_prob_results/third_order_mu0.01/fixation_probs_15.dat')
else:
    # we are in a python script
    # input file is specified as a command line argument
    if len(sys.argv) < 2:
        print("Usage: python trans_graph_vis.py <input_file>")
        sys.exit(1)
    input_file = sys.argv[1]

input_file
# %%
# Open the file for reading

norm_ids = []
self_coop_level = []
p_fix = []

with open(input_file, 'r') as file:
    # Iterate through each line in the file
    for line in file:
        if line.startswith('#'):
            continue
        # Split the line into individual values
        values = line.split()

        # Parse the first two values and append them to the respective lists
        norm_ids.append(int(values[0]))
        self_coop_level.append(float(values[1]))

        # Parse the remaining values and append them to the data_matrix list
        p_fix.append([float(val) for val in values[2:]])

norm_ids = np.array(norm_ids)
self_coop_level = np.array(self_coop_level)
p_fix = np.array(p_fix)

norm_ids, self_coop_level, p_fix

# %%
def extract_strategies(norm_ids, self_coop_level, p_fix, strategies):
    idx = [np.where(norm_ids==s)[0][0] for s in strategies]
    return norm_ids[idx], self_coop_level[idx], p_fix[idx][:,idx]

names = ["ALLD", "L1", "L2", "L3", "L4", "L5", "L7"]
selection = [64704, 765131, 634059, 769226, 761034, 638154, 859333]
selected_norm_ids, selected_pc, selected_p_fix = extract_strategies(norm_ids, self_coop_level, p_fix, selection)

selected_norm_ids, selected_pc, selected_p_fix

# %%
import networkx as nx
from pyvis.network import Network

G = nx.DiGraph()

n = len(selected_norm_ids)
for i in range(n):
    G.add_node(i, norm_id=names[i])
for i in range(n):
    for j in range(n):
        if i == j:
            continue
        if selected_p_fix[i,j] > 0.02/4:
            G.add_edge(i, j, weight=selected_p_fix[i,j])

G

# %%
plt.clf()
pos = nx.spring_layout(G)  # You can choose a different layout if you prefer
edge_labels = {(i, j): G[i][j]['weight'] for i, j in G.edges()}
node_labels = {i: G.nodes[i]['norm_id'] for i in G.nodes()}
nx.draw(G, pos, with_labels=False, node_color='skyblue', node_size=1500, font_size=12, font_weight='bold', width=2)
nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=10)
nx.draw_networkx_labels(G, pos, labels=node_labels, font_size=10)

# Show the graph
plt.show()

# %%
# visualize G
nt = Network(notebook=True, cdn_resources='remote', directed=True, select_menu=True)
for node, norm_id in nx.get_node_attributes(G, 'norm_id').items():
    nt.add_node(node, label=str(norm_id))
    print(norm_id)
for edge in G.edges(data=True):
    nt.add_edge(edge[0], edge[1], value=edge[2]['weight'])  # Include edge weights
    #label = '%.2f' % (edge[2]['weight'],)
    #nt.add_edge(edge[0], edge[1], value=edge[2]['weight'], label=label)  # Include edge weights

nt.toggle_physics(False)

#nt.show_buttons(filter_=['physics', 'nodes', 'edges', 'layout', 'interaction', 'manipulation', 'selection', 'renderer'])
nt.show_buttons(filter_=['physics'])
nt.show('nx.html', notebook=True)
#print(nt.template)
# %%

def show_norms(norm_ids_s):
    norm_names = {
        # 64704: "ALLD",
        765131: "L1",
        634059: "L2",
        769226: "L3",
        761034: "L4",
        638154: "L5",
        629962: "L6",
        859333: "L7",  # or 756938
        892101: "L8", # or 625866
        765130: "L1 BBD",
        765129: "L1 BGD",
        634058: "L2 BBD",
        634057: "L2 BGD",
        634050: "L2 GGD BBD",
        769227: "L3 BBC",
        761035: "L4 BBC",
        638155: "L5 BBC",
        859341: "L7 BBC",
        568523: "L2 GBDB"
    }
# %%
