import networkx as nx
import pandas as pd
import numpy as np
import sys, networkx as nx, matplotlib.pyplot as plt
import argparse
from math import sqrt
import math
import statistics
import datetime
import csv
import pandas as pd
import numpy as np
import re
from collections import Counter
import glob
import os
import scipy.stats as stats
from scipy.stats import ttest_ind
### TO remove future warnings
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

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
    
    parser.add_argument("general_path")
    parser.add_argument("disease_path")
    parser.add_argument("VDJdb_path")
    parser.add_argument("output_path")

    parser.add_argument('--sourceSequence')
    options = parser.parse_args()
    return options
# setting the path for joining multiple files

if __name__ == "__main__":

	options = getArguments()
	file_name = options.general_path
	diseasePath = options.disease_path
	VDJPath = options.VDJdb_path
	plotsPath = options.output_path
	source = options.sourceSequence

def Average(lst):
    return int(sum(lst) / len(lst))

if 'TRA' in file_name:
    chain = 'TRA'
elif 'TRB' in file_name:
    chain = 'TRB'

print(chain)

df = pd.read_csv(file_name, low_memory = False)
a = pd.read_csv(diseasePath, low_memory = False)
b = pd.read_csv(VDJPath,low_memory=False)
plots_saving_path = plotsPath




junctions = []

d = pd.DataFrame()

junctions = a['Junction']

junctions = list(set(junctions))

qq = pd.DataFrame()


for junction in junctions:
    if junction == 'Junction' or pd.isnull(junction):
        continue
    else:
        q = a[a['Junction']==junction]
        qq = qq.append(q)
        
        c = b[b['CDR3']==junction]
        if not c.empty:
            

            d = d.append(c,ignore_index=True)
            d = d[['CDR3','V','J','Species','Epitope','Epitope species']]

            graph = pd.DataFrame()
            graph['diseases'] = (a[a['Junction']==junction]['Disease'])
            graph['junction_hits'] = (a[a['Junction']==junction]['Total_hits'])
            graph['junction_hits'] = graph['junction_hits'].astype('int')
            
            graph = graph.sort_values('junction_hits').reset_index(drop=True)
            
            plt.barh(graph['diseases'], graph['junction_hits'])
            
            plt.xlabel('junction_hits')
            plt.title("{0} ----> {1}(Epitope) ----> {2} : Found in VDJdb".format(junction,c['Epitope'].values[0],c['Epitope species'].values[0]))
            plt.tight_layout()
            plt.savefig('{0}{1}_{2}'.format(plots_saving_path,junction,chain),dpi=200, bbox_inches = "tight")
            plt.show()
            
            
            
            print('\nFOUNDDD FOR {0}\n'.format(junction))
            

qq['Total_hits'] = qq['Total_hits'].astype('int')


if d.empty:
    print('No matches found in VDJdb')
else:

    print(d)
    #stat, p_value = ttest_ind(graph['junction_hits'], qq['Total_hits'])
    #print(f"\nt-test: statistic={stat:.4f}, p-value={p_value:.4f}\n")


distribution_all = {}


#for d in set(qq['Disease']):
    #hits = []
    #for j in qq['Junction']:
        #hit = qq[qq['Disease']==d]['Total_hits'].values[0]
        #print(hit)
        #hits.append(hit)
        
    #distribution_all[d] = hits

def Average(lst):
    return sum(lst) / len(lst)     

for d in set(qq['Disease']):
    hits = []
    for j in set(qq['Junction']):
        subset = qq[(qq['Disease'] == d) & (qq['Junction'] == j)]
        if not subset.empty:
            
            hit = subset['Total_hits'].values[0]
            hits.append(hit)
            
    distribution_all[d] = hits  

important_sequences = []

for sequence in set(qq['Junction']):
    if df[df['Sequence']== sequence]['All_hits'].values[0] > Average(df['All_hits']):
        important_sequences.append(sequence)
    else:
        continue
important_sequences_df = pd.DataFrame()

for j in important_sequences:
    important = qq[qq['Junction'] == j]
    important_sequences_df = important_sequences_df.append(important)
    

distribution_plot = {k: sum(distribution_all[k]) for k in distribution_all.keys()}

x_pie = [t for t in distribution_plot.values()]
labels_pie = [z for z in distribution_plot.keys()]
fig, ax = plt.subplots()
wedges, texts, autotexts = ax.pie(x_pie,labels = labels_pie,autopct='%1.1f%%' , startangle = 90)
ax.legend(bbox_to_anchor=(1.0005, 1.0), loc='upper left')

ax.set_title('All hits')


threshold = 10
for label, pct_label in zip(texts, autotexts):
    pct_value = pct_label.get_text().rstrip('%')
    pct_value = int(float(pct_value))
    print(pct_value)
    if float(pct_value) < threshold:
        label.set_text('')
        pct_label.set_text('')

plt.savefig('{0}{1}_{2}'.format(plots_saving_path,'PieChart',chain),dpi=300, bbox_inches = "tight")
plt.show()


# plt.pie(distribution_plot.values(), labels=distribution_plot.keys(),
#         startangle=90, autopct=lambda p: format(p, '.2f') if p > 25 else None, colors=plt.cm.Set2.colors)
# plt.title('All hits')
# plt.tight_layout()
# plt.show()

print(len(set(qq['Disease'])))

#print("\nHERREE{0}".format(distribution_all['breast cancer']))

for i in set(qq['Disease']):
    
    plt.title('{0}'.format(i))
    hist, bins = np.histogram(distribution_all[i], bins=50)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    plt.bar(center, hist, align='center', width=width)
    plt.axvline(Average(distribution_all[i]), color='red', linewidth=1)
    plt.xlabel('# of patients', fontsize = 10)
    plt.ylabel('# of sequences', fontsize = 10)

    q90 = np.quantile(distribution_all[i], 0.9)
    plt.axvline(q90, color='yellow', linewidth=1)

    plot_text = {}

    for j in (important_sequences_df[important_sequences_df['Disease'] == i]['Total_hits']):
        plot_text[j] = ''
    
    for index, j in enumerate(important_sequences_df[important_sequences_df['Disease'] == i]['Total_hits']):
           
        #plt.axvline(j, color='green', linestyle='dashed', linewidth=1)
        
        #print("FIRST{0}".format(j - Average(distribution_all[i])))
        
        #print("SECOND{0}".format(max ((important_sequences_df[important_sequences_df['Disease'] == i]['Total_hits']) - Average(distribution_all[i]))))
        
        if j > q90:
            sequencesss = ("'{0}' ".format(important_sequences_df[important_sequences_df['Disease'] == i]['Junction'].iloc[index]))
            #print(sequencesss)
            plot_text[j]+=(sequencesss)
            
            print
           
            #plot_text[j]=(important_sequences_df[(important_sequences_df['Disease'] == i) & (important_sequences_df['Total_hits'] == j)]['Junction'].values[0])
    for j in plot_text:
        print(plot_text[j])
        plt.text(j, 0.9,plot_text[j],rotation = 270,size = 5,fontstretch = 'extra-expanded',color= 'red')
        
    #plt.tight_layout()
    plt.savefig('{0}{1}_{2}'.format(plots_saving_path,i,chain),dpi=150, bbox_inches = "tight")

    plt.show()























nodes = df['Sequence']
nodes = list(nodes)


node_sizes = df['All_hits']
node_sizes = list(node_sizes)
print(df[df['Sequence']==source]['All_hits'])

print(len(nodes))
print(len(node_sizes))
labels = {}
for n in nodes:
    if df[df['Sequence']== n]['All_hits'].values[0] > Average(node_sizes)*1.2 :
    	labels[n] =  '{0} : {1} CD4; {2} CD8 ; {3} patient'.format(n,(df[df['Sequence']== n]['CD4_count'].values[0]), df[df['Sequence']== n]['CD8_count'].values[0],df[df['Sequence']== n]['Patients'].values[0])

    elif df[df['Sequence']== n]['All_hits'].values[0] > Average(node_sizes):
    	labels[n] = n

    else:
    	labels[n] = ''

# Node sizes: [0, 100, 200, 300, 400, 500, 600, 700, 800, 900]

# Connect each node to its successor
edges = [ ('{0}'.format(source), i) for i in nodes if i != '{0}'.format(source) ]

# Create the graph and draw it with the node labels
g = nx.Graph()
g.add_nodes_from(nodes)
g.add_edges_from(edges)
color_map = []
for node in node_sizes:
    if node > statistics.median(node_sizes):
        color_map.append('red')
    else:
        color_map.append('blue')


plt.suptitle('{0} (Wild-type sequence) has {1} hits and {2} different mutated sequences'.format(source,df[df['Sequence']==source]['All_hits'].values[0],len(node_sizes)-1))

plt.title('All hits: {0}    Total CD4 count: {1}    Total CD8 count: {2}'.format(df['All_hits'].sum(),df['CD4_count'].sum(),df['CD8_count'].sum()), y=-0.1)

nx.draw_networkx(g, node_size = node_sizes,node_color='red',edge_color= 'black',width = 0.1, labels=labels,font_size = 7, with_labels=True,alpha = 0.6,font_weight = 'normal')
for n in [min(node_sizes),math.trunc(Average(node_sizes)),max(node_sizes)]:
    plt.plot([], [], 'bo', markersize = sqrt(n), label = f"{n}",color= 'black',alpha = 0.4)
leg = plt.legend(fontsize=13,labelspacing = 5, loc='center left', bbox_to_anchor=(1, 0.5), frameon = False)
for text in leg.get_texts():
    text.set_color("black")

plt.box(False)
#plt.savefig('{0}All_sequences_visualization'.format(plots_saving_path),dpi=450, bbox_inches = "tight")

plt.show()

