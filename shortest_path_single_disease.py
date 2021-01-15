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

url_to_common_name = {line.strip().split()[0]:line.strip().split()[1] for line in open('Edgelists/PKT_Master_Edge_List_NoOntologyData_LABELS.txt')}

G = nx.Graph()
# for each disease gene set
for line in open('Edgelists/PKT_Master_Edge_List_NoOntologyData.txt'):
    edge = line.strip().split('\t')
    if len(edge) != 4:
        continue
    subject = edge[1]
    the_object = edge[3]
    if subject in url_to_common_name:
        subject = url_to_common_name[subject]
    if the_object in url_to_common_name:
        the_object = url_to_common_name[the_object]
    G.add_edge(subject, the_object, predicate=edge[2], predicate_name=edge[0])

gene_sets = 'sIMB_and_random_shuffle_disease_set.tsv'
disease_gene_sets = {}
for line in open(gene_sets, 'r'):
    row = line.strip().split('\t')
    # the first item is the disease common name, everything else is genes
    disease_gene_sets[row[0]] = row[1:]


if True or not os.path.exists('sIMB_node_counts.pickle') or not os.path.exists('sIMB_counts.pickle') or not os.path.exists('sIMB_path_lengths.pickle'):
    counts = {}
    node_counts = {}
    path_lengths = {'length': [], 'type': []}
    d_keys = list(disease_gene_sets.keys())

    disease = 'Inclusion Body Myopathy, Sporadic'
    disease_genes = disease_gene_sets[disease]
    nodes = partition(disease_genes, 2)
    seeds, targets = nodes[0], nodes[1]
    for i, seed in enumerate(seeds):
        for target in targets:
            try:
                path = nx.shortest_path(G, seed, target)
            except (nx.exception.NodeNotFound, nx.exception.NetworkXNoPath):
                continue
            path_lengths['length'].append(len(path))
            path_lengths['type'].append('normal')
            for i in range(len(path) - 1):
                p = G[path[i]][path[i + 1]]['predicate_name']
                # add proportional weight
                if p in counts:
                    counts[p] += 1
                else:
                    counts[p] = 1
            # count the number times nodes occur in the paths
            for node in path:
                if node in node_counts:
                    node_counts[node] += 1
                else:
                    node_counts[node] = 1
    pickle.dump(counts, open('sIMB_counts.pickle', 'wb'))
    pickle.dump(node_counts, open('sIMB_node_counts.pickle', 'wb'))
    pickle.dump(path_lengths, open('sIMB_path_lengths.pickle', 'wb'))

counts = pickle.load(open('sIMB_counts.pickle', 'rb'))
node_counts = pickle.load(open('sIMB_node_counts.pickle', 'rb'))
path_lengths = pickle.load(open('sIMB_path_lengths.pickle', 'rb'))

df = pd.DataFrame({'predicate':list(counts.keys()),'count':list(counts.values())})
df = df.sort_values('predicate')
g = sns.barplot(df['predicate'], df['count'])
g.set_yscale("log")
plt.xticks(rotation=90)
plt.tight_layout()
plt.ylabel('Count')
plt.xlabel('Predicate')
plt.savefig('Figures/sIMB_predicate_counts_log.png')
plt.show()

g = sns.barplot(df['predicate'], df['count'])
plt.xticks(rotation=90)
plt.tight_layout()
plt.ylabel('Count')
plt.xlabel('Predicate')
plt.savefig('Figures/sIMB_predicate_counts.png')
plt.show()

# do it again with shuffled targets


if True or not os.path.exists('shuffled_sIMB_node_counts.pickle') or not os.path.exists('shuffled_sIMB_counts.pickle') or not os.path.exists('shuffled_sIMB_path_lengths.pickle'):
    shuffled_counts = {}
    shuffled_node_counts = {}
    shuffled_path_lengths = {'length': [], 'type': []}
    d_keys = list(disease_gene_sets.keys())

    disease = 'Inclusion Body Myopathy, Sporadic'
    disease_genes = disease_gene_sets[disease]
    nodes = partition(disease_genes, 2)
    seeds, targets = nodes[0], nodes[1]
    shuffled_nodes = partition(disease_gene_sets['Randomly Chosen Null Set'], 2)
    shuffled_seeds, shuffled_targets = shuffled_nodes[0], shuffled_nodes[1]
    for i, seed in enumerate(seeds):
        # use shuffled targets
        for target in shuffled_targets:
            try:
                path = nx.shortest_path(G, seed, target)
            except (nx.exception.NodeNotFound, nx.exception.NetworkXNoPath):
                continue

            shuffled_path_lengths['length'].append(len(path))
            shuffled_path_lengths['type'].append('shuffled')
            for i in range(len(path) - 1):
                p = G[path[i]][path[i + 1]]['predicate_name']
                # add proportional weight
                if p in shuffled_counts:
                    shuffled_counts[p] += 1
                else:
                    shuffled_counts[p] = 1
            # count the number times nodes occur in the paths
            for node in path:
                if node in shuffled_node_counts:
                    shuffled_node_counts[node] += 1
                else:
                    shuffled_node_counts[node] = 1
    pickle.dump(shuffled_counts, open('shuffled_sIMB_counts.pickle', 'wb'))
    pickle.dump(shuffled_node_counts, open('shuffled_sIMB_node_counts.pickle', 'wb'))
    pickle.dump(shuffled_path_lengths, open('shuffled_sIMB_path_lengths.pickle', 'wb'))

shuffled_counts = pickle.load(open('shuffled_sIMB_counts.pickle', 'rb'))
shuffled_node_counts = pickle.load(open('shuffled_sIMB_node_counts.pickle', 'rb'))
shuffled_path_lengths = pickle.load(open('shuffled_sIMB_path_lengths.pickle', 'rb'))

shuffled_df = pd.DataFrame({'predicate':list(shuffled_counts.keys()),'count':list(shuffled_counts.values())})
shuffled_df = shuffled_df.sort_values('predicate')
g = sns.barplot(shuffled_df['predicate'], df['count'])
g.set_yscale("log")
plt.xticks(rotation=90)
plt.tight_layout()
plt.ylabel('Count')
plt.xlabel('Predicate')
plt.savefig('Figures/shuffled_sIMB_predicate_counts_log.png')
plt.show()

g = sns.barplot(shuffled_df['predicate'], df['count'])
plt.xticks(rotation=90)
plt.tight_layout()
plt.ylabel('Count')
plt.xlabel('Predicate')
plt.savefig('Figures/shuffled_sIMB_predicate_counts.png')
plt.show()

# plot change
# get all predicates
all_preds = set(df['predicate']).union(set(shuffled_df['predicate']))
diff_dict = {'predicate': [], 'difference': []}
for p in all_preds:
    norm_count = 0
    shuf_count = 0
    if p in counts:
        norm_count = counts[p]
    if p in shuffled_counts:
        shuf_count = shuffled_counts[p]
    diff_dict['predicate'].append(p)
    diff_dict['difference'].append(norm_count - shuf_count)

sns.barplot(diff_dict['predicate'],diff_dict['difference'])
plt.xticks(rotation=90)
plt.tight_layout()
plt.ylabel('Difference (Normal Count - Shuffled Count)')
plt.xlabel('Predicate')
plt.savefig('Figures/sIMB_norm_shuffled_change.png')
plt.show()

# plot path lengths:
length_df = pd.concat([pd.DataFrame(path_lengths),pd.DataFrame(shuffled_path_lengths)])
sns.boxplot(length_df['type'],length_df['length'])
plt.savefig('Figures/sIMB_norm_shuffled_length_boxplots.png')
plt.show()

filtered_node_counts = {n:node_counts[n] for n in node_counts if n not in disease_gene_sets['Inclusion Body Myopathy, Sporadic']}
node_counts_df = pd.DataFrame({'node':filtered_node_counts.keys(),'count':filtered_node_counts.values()}).sort_values('count',ascending=False)

filtered_shuffled_node_counts = {n:shuffled_node_counts[n] for n in shuffled_node_counts if n not in disease_gene_sets['Randomly Chosen Null Set']}
shuffled_node_counts_df = pd.DataFrame({'node':filtered_shuffled_node_counts.keys(),'count':filtered_shuffled_node_counts.values()}).sort_values('count',ascending=False)
