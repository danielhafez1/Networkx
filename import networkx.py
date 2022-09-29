import networkx as nx
import pandas as pd
import numpy as np
import sys, networkx as nx, matplotlib.pyplot as plt
import argparse

def getArguments():
    # Set up the command line parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="")

    parser.add_argument(
    "-v",
    "--verbose",
    action="store_true",
    help="Run the program in verbose mode.")
    
    parser.add_argument("path")
    parser.add_argument('--sourceSequence')
    options = parser.parse_args()
    return options
# setting the path for joining multiple files

if __name__ == "__main__":

	options = getArguments()
	file_name = options.path
	source = options.sourceSequence




df = pd.read_csv(file_name, low_memory = False)

nodes = df['Sequence']
nodes = list(nodes)


node_sizes = df['All_hits']
node_sizes = list(node_sizes)


print(len(nodes))
print(len(node_sizes))
labels = {}
for n in nodes:
    
    labels[n] =  n

# Node sizes: [0, 100, 200, 300, 400, 500, 600, 700, 800, 900]

# Connect each node to its successor
edges = [ ('{0}'.format(source), i) for i in nodes if i != '{0}'.format(source) ]

# Create the graph and draw it with the node labels
g = nx.Graph()
g.add_nodes_from(nodes)
g.add_edges_from(edges)

try:

	nx.draw(g, node_size = node_sizes, labels=labels,font_size = 7, with_labels=True,alpha = 0.5,font_weight = 'normal')    
	plt.show()

except:
	print("Network draw doesnt work due to the absence of the right source sequence")












