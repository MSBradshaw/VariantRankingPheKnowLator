import networkx as nx
import pickle
import random
import rdflib
import re
import pandas as pd
from rwr import partition
import os
import matplotlib.pyplot as plt
import seaborn as sns

G = nx.Graph()
# for each disease gene set
for line in open('Edgelists/pkl_edgelist_gene_symbols_no_snps.tsv'):
    edge = line.strip().split('\t')
    if len(edge) != 3:
        continue
    G.add_edge(edge[0], edge[2], predicate=edge[1])

gene_sets = 'DisGeNET_genesets.txt'
disease_gene_sets = {}
for line in open(gene_sets, 'r'):
    row = line.strip().split('\t')
    # the first item is the disease common name, everything else is genes
    disease_gene_sets[row[0]] = row[1:]

if not os.path.exists('node_counts.pickle') or not os.path.exists('counts.pickle') or not os.path.exists('path_lengths.pickle'):
    counts = {}
    node_counts = {}
    path_lengths = {'length': [], 'type': []}
    d_keys = list(disease_gene_sets.keys())
    for i, disease in enumerate(disease_gene_sets.keys()):
        print(i)
        disease_genes = disease_gene_sets[disease]
        nodes = partition(disease_genes, 2)
        seeds, targets = nodes[0], nodes[1]
        for i, seed in enumerate(seeds):
            for target in targets:
                try:
                    path = nx.shortest_path(G, seed, target)
                except nx.exception.NodeNotFound:
                    pass
                if 'http://purl.obolibrary.org/obo/SO_0001217' in path or 'http://purl.obolibrary.org/obo/NCBITaxon_9606' in path:
                    path_lengths['length'].append(len(path))
                    path_lengths['type'].append('outlier')
                else:
                    path_lengths['length'].append(len(path))
                    path_lengths['type'].append('normal')
                for i in range(len(path) - 1):
                    p = G[path[i]][path[i + 1]]['predicate']
                    # add proportional weight
                    if p in counts:
                        counts[p] += 1 / len(path)
                    else:
                        counts[p] = 1 / len(path)
                # count the number times nodes occur in the paths
                for node in path:
                    if node in node_counts:
                        node_counts[node] += 1
                    else:
                        node_counts[node] = 1
    pickle.dump(counts, open('counts.pickle', 'wb'))
    pickle.dump(node_counts, open('node_counts.pickle', 'wb'))
    pickle.dump(path_lengths, open('path_lengths.pickle', 'wb'))

counts = pickle.load(open('counts.pickle', 'rb'))
node_counts = pickle.load(open('node_counts.pickle', 'rb'))
path_lengths = pickle.load(open('path_lengths.pickle', 'rb'))

# plot the counts (they will be proportional)
# are there some over represented nodes?
# do these over represented nodes have over represented edge types
plt.hist(list(node_counts.values()))
plt.yscale('log')
plt.xlabel('Occurrences')
plt.ylabel('Node Count')
plt.savefig('Figures/node_occurrences.png')
plt.show()

df = pd.DataFrame({'node': list(node_counts.keys()), 'count': list(node_counts.values())}).sort_values('count',
                                                                                                       ascending=False)
print(df)

# get all the edges types connected to
list(G.neighbors('http://purl.obolibrary.org/obo/SO_0001217'))

predicates = {}
for neighbor in G.neighbors('http://purl.obolibrary.org/obo/SO_0001217'):
    pred = G.edges['http://purl.obolibrary.org/obo/SO_0001217', neighbor]['predicate']
    if pred in predicates:
        predicates[pred] += 1
    else:
        predicates[pred] = 1
print(predicates)

predicates = {}
for neighbor in G.neighbors('http://purl.obolibrary.org/obo/NCBITaxon_9606'):
    pred = G.edges['http://purl.obolibrary.org/obo/NCBITaxon_9606', neighbor]['predicate']
    if pred in predicates:
        predicates[pred] += 1
    else:
        predicates[pred] = 1
print(predicates)

name_mapping = {'http://purl.obolibrary.org/obo/RO_0000056': 'participates_in',
                'http://purl.obolibrary.org/obo/RO_0002160': 'evolutionarily related to',
                'http://www.w3.org/2000/01/rdf-schema#subClassOf': 'sub class of',
                'http://purl.obolibrary.org/obo/RO_0003302': 'causes or contributes to condition',
                'http://purl.obolibrary.org/obo/RO_0002511': 'Inverse of transcribed from',
                'http://purl.obolibrary.org/obo/RO_0001025': 'located in',
                'http://purl.obolibrary.org/obo/RO_0002205': 'has gene product',
                'http://purl.obolibrary.org/obo/RO_0002606': 'is substance that treats',
                'http://purl.obolibrary.org/obo/RO_0002434': 'interacts with',
                'http://purl.obolibrary.org/obo/RO_0002200': 'has phenotype',
                'http://purl.obolibrary.org/obo/BFO_0000050': 'part of',
                'http://purl.obolibrary.org/obo/RO_0000052': 'inheres in',
                'http://purl.obolibrary.org/obo/CLO_0000179': ' is disease model for',
                'http://purl.obolibrary.org/obo/RO_0002435': 'genetically interacts with',
                'http://purl.obolibrary.org/obo/RO_0004019': 'disease has basis in',
                'http://purl.obolibrary.org/obo/CLO_0000015': 'derives from patient having disease',
                'http://purl.obolibrary.org/obo/RO_0002436': 'molecularly interacts with',
                'http://purl.obolibrary.org/obo/RO_0002573': 'has modifier',
                'http://purl.obolibrary.org/obo/IDO_0000664': 'has_material_basis_in',
                'http://purl.obolibrary.org/obo/so#has_part': 'has part',
                'http://purl.obolibrary.org/obo/so#non_functional_homolog_of': 'non_functional_homolog_of',
                'http://purl.obolibrary.org/obo/so#member_of': 'member of',
                'http://purl.obolibrary.org/obo/RO_0002314': 'inheres in part of',
                'http://purl.obolibrary.org/obo/CLO_0000167': 'has disease',
                'http://purl.obolibrary.org/obo/RO_0000086': 'has quality',
                'http://purl.obolibrary.org/obo/so#part_of': 'part of',
                'http://purl.obolibrary.org/obo/BFO_0000051': 'has part'}

named_counts = {name_mapping[x]: counts[x] for x in counts.keys()}

plotting_df = pd.DataFrame({'name': list(named_counts.keys()), 'count': list(named_counts.values())})
plotting_df = plotting_df.sort_values('name')
# plt.hist(df['name'], df['count'])


g = sns.barplot(plotting_df['name'], plotting_df['count'])
g.set_yscale("log")
plt.xticks(rotation=90)
plt.tight_layout()
plt.ylabel('Proportional Count')
plt.xlabel('Predicate')
plt.savefig('Figures/proportional_predicate_counts.png')
plt.show()

# box plots of path lengths among outliers and none outliers using path_lengths dictionary
length_df = pd.DataFrame(path_lengths)

sns.boxplot('type','length', data=length_df)
plt.ylabel('Path Length')
plt.xlabel(' ')
plt.savefig('Figures/path_length.png')
plt.show()
