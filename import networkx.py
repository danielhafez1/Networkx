import networkx as nx
import pandas as pd
import numpy as np
import sys, networkx as nx, matplotlib.pyplot as plt


filepath = '/Users/daniel/Desktop/Master-Project/ireceptor_script-main/Camilla_results/General_SameVJ_I-receptor_TRB.csv'

df = pd.read_csv(filepath, low_memory = False)

nodes = df['Sequence']
nodes = list(nodes)
print(nodes)
print(len(nodes))

node_sizes = df['All_hits']
node_sizes = list(node_sizes)
print(node_sizes)
print(len(node_sizes))



labels = {}
for n in nodes:
    
    labels[n] =  n

# Node sizes: [0, 100, 200, 300, 400, 500, 600, 700, 800, 900]

# Connect each node to its successor
edges = [ ('CSARRGGNTGELFF', i) for i in nodes if i != 'CSARRGGNTGELFF' ]

# Create the graph and draw it with the node labels
g = nx.Graph()
g.add_nodes_from(nodes)
g.add_edges_from(edges)

nx.draw(g, node_size = node_sizes, labels=labels,font_size = 7, with_labels=True,alpha = 0.5,font_weight = 'normal')    
plt.show()












