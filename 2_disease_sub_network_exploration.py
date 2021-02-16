import networkx as nx
from rwr import partition
from webweb import Web
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


"""
--edge_list
Edgelists/PKT_Master_Edge_List_NoOntologyData.txt
--disease_sets
sIMB_only_disease_set.txt
--node_mapping
Edgelists/PKT_Master_Edge_List_NoOntologyData_LABELS.txt
--threshold
0.0000165
"""
node_mapping = 'Edgelists/PKT_Master_Edge_List_NoOntologyData_LABELS.txt'
el = 'Edgelists/PKT_Master_Edge_List_NoOntologyData.txt'

url_to_common_name = {line.strip().split()[0]: line.strip().split()[1] for line in open(node_mapping)}

G = nx.Graph()
types = set()
for line in open(el, 'r'):
    row = line.strip().split('\t')
    if len(row) == 3:
        # when using 'Edgelists/pkl_triples_with_symbols.tsv'
        G.add_edge(row[0], row[2])
        G.edges[row[0], row[2]]['type'] = row[1]
        G.edges[row[0], row[2]]['weight'] = .001
        types.add(row[0])
    elif len(row) == 4:
        subject = row[1]
        the_object = row[3]
        if subject in url_to_common_name:
            subject = url_to_common_name[subject]
        if the_object in url_to_common_name:
            the_object = url_to_common_name[the_object]
        # when using 'Edgelists/PKT_Master_Edge_List_NoOntologyData.txt'
        G.add_edge(subject, the_object)
        G.edges[subject, the_object]['type'] = row[2]
        G.edges[subject, the_object]['type_name'] = row[0]
        G.edges[subject, the_object]['weight'] = .001
        types.add(row[0])

# get all shortest paths from both diseases

gene_sets = '2_disease.txt'
disease_gene_sets = {}
for line in open(gene_sets, 'r'):
    row = line.strip().split('\t')
    # the first item is the disease common name, everything else is genes
    disease_gene_sets[row[0]] = row[1:]

paths = {}
for disease in disease_gene_sets.keys():
    print(disease)
    paths[disease] = set()
    disease_genes = disease_gene_sets[disease]
    nodes = partition(disease_genes, 2)
    seeds, targets = nodes[0], nodes[1]
    for i, seed in enumerate(seeds):
        for target in targets:
            try:
                temp_paths = nx.all_shortest_paths(G, seed, target)
            except (nx.exception.NodeNotFound, nx.exception.NetworkXNoPath):
                continue
            for path in temp_paths:
                paths[disease] = paths[disease].union(set(path))

# get a sub network from each of these genes
subnetworks = {}
for disease in paths.keys():
    nodes = list(paths[disease])
    g = G.subgraph(nodes)
    subnetworks[disease] = g


for disease in subnetworks.keys():
    display = {
        'nodes': {str(i): {'name': node,'deg':G.degree(node), 'is seed or target': node in disease_gene_sets[disease]} for i, node in enumerate(subnetworks[disease])}
    }
    # web = Web(adjacency_matrix, adjacency_type='matrix')
    web = Web(nx.to_numpy_array(subnetworks[disease]),display=display)
    # show the visualization
    web.display.colorBy = 'is seed or target'
    web.display.sizeBy = 'deg'
    web.display.showNodeNames = True
    web.display.linkLength = 100
    web.display.charge = 300
    web.display.radius = 10
    web.show()


# --------------------------------------- Huntington's Disease: HTT ---------------------------------------
# all neighbors of HTT
len(list(G.neighbors('HTT')))
# get the predicate count of edges articulating CFTR
htt_predicates = {}
for node in G.neighbors('HTT'):
    edge = G.edges[node,'HTT']
    if edge['type_name'] in htt_predicates:
        htt_predicates[edge['type_name']] += 1
    else:
        htt_predicates[edge['type_name']] = 1

htt_df = pd.DataFrame({'predicate':list(htt_predicates.keys()),'count':list(htt_predicates.values())}).sort_values('count',ascending=False)

g = sns.barplot(htt_df['predicate'], htt_df['count'])
plt.xticks(rotation=90)
plt.tight_layout()
plt.ylabel('Count')
plt.xlabel('Predicate')
plt.savefig('Figures/HTT_predicate_counts.png')
plt.show()

# --------------------------------------- Cystic Fibrosis: HTT ---------------------------------------
len(list(G.neighbors('CFTR')))
# get the predicate count of edges articulating CFTR
cftr_predicates = {}
for node in G.neighbors('CFTR'):
    edge = G.edges[node,'CFTR']
    if edge['type_name'] in htt_predicates:
        cftr_predicates[edge['type_name']] += 1
    else:
        cftr_predicates[edge['type_name']] = 1

cftr_df = pd.DataFrame({'predicate':list(cftr_predicates.keys()),'count':list(cftr_predicates.values())}).sort_values('count',ascending=False)

g = sns.barplot(cftr_df['predicate'], cftr_df['count'])
plt.xticks(rotation=90)
plt.tight_layout()
plt.ylabel('Count')
plt.xlabel('Predicate')
plt.savefig('Figures/CFTR_predicate_counts.png')
plt.show()

# get diameter of network
# remove non-full subgraphs
connected_nodes = None
for x in nx.connected_components(G):
    if len(x) > 500:
        print('Greater')
        connected_nodes = x

nx.diameter(nx.subgraph(G,connected_nodes))