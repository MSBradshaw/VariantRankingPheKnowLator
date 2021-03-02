import networkx as nx
from webweb import Web
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import obonet
from rwr import partition
import numpy as np
import pickle
import os

url = 'https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.obo'
hpo_g = obonet.read_obo(url)

G = nx.read_edgelist('HPO_String_Analysis/HPO_String_edgelist.tsv', delimiter='\t')
# HPO = nx.read_edgelist('HPO_String_Analysis/HPO_String_edgelist.tsv', delimiter='\t')

# remove connection between CFTR and it's known HPO terms
terms_to_remove = ['HP:0011949', 'HP:0001733', 'HP:0002202', 'HP:0030830', 'HP:0000006', 'HP:0011962', 'HP:0002105',
                   'HP:0000007', 'HP:0004313', 'HP:0002027', 'HP:0002721', 'HP:0005952', 'HP:0006528', 'HP:0006538',
                   'HP:0002110', 'HP:0030877', 'HP:0002613', 'HP:0011961', 'HP:0001738', 'HP:0030247', 'HP:0002783',
                   'HP:0000118', 'HP:0000837', 'HP:0005376', 'HP:0001508', 'HP:0001658', 'HP:0003251', 'HP:0000819',
                   'HP:0012873', 'HP:0100749', 'HP:0002099', 'HP:0005213', 'HP:0002024', 'HP:0004401', 'HP:0000798',
                   'HP:0001425', 'HP:0002205', 'HP:0100027', 'HP:0012379', 'HP:0001945', 'HP:0004469', 'HP:0011947',
                   'HP:0002240', 'HP:0002150', 'HP:0011227', 'HP:0001974', 'HP:0004326', 'HP:0002570', 'HP:0001648',
                   'HP:0002035']

cftr_hpos = terms_to_remove

for hpo in terms_to_remove:
    try:
        G.remove_edge(hpo, 'CFTR')
    except nx.exception.NetworkXError:
        print('error')
        pass

# -------------------------------------

paths = []
# length of paths getting from those terms to CFTR now
for hpo in terms_to_remove:
    try:
        temp_paths = list(nx.all_shortest_paths(G, hpo, 'CFTR'))
    except:
        print('No path from ' + hpo + ' to CFTR')
        continue
    if len(list(temp_paths)[0]) > 3:
        print(hpo)
    # print(temp_paths)
    for path in temp_paths:
        # print(len(path))
        # show the visualization
        paths.append(path)

nodes = []
node_counts = {}
for path in paths:
    for p in path:
        nodes.append(p)
        if p in node_counts:
            node_counts[p] += 1
        else:
            node_counts[p] = 1

g = nx.subgraph(G, nodes)

display = {
    'nodes': {str(i): {'name': node, 'deg': G.degree(node), 'related': node in terms_to_remove,
                       'node_count': node_counts[node]} for i, node in enumerate(g.nodes())}
}

web = Web(nx.to_numpy_array(g), display=display)
web.display.colorBy = 'related'
web.display.sizeBy = 'deg'
web.display.showNodeNames = True
web.display.linkLength = 100
web.display.charge = 300
web.display.radius = 10
web.show()

keys = list(node_counts.keys())
filtered_counts = [node_counts[key] for key in keys if node_counts[key] > 50 and key != 'CFTR']
filtered_keys = [key for key in keys if node_counts[key] > 50 and key != 'CFTR']

plt.hist(filtered_counts)
plt.show()

named_keys = ['Autosomal dominant inheritance', 'Autosomal recessive inheritance', 'Abdominal pain',
              'Phenotypic abnormality', 'Abnormality of metabolism/homeostasis', 'Failure to thrive', 'Heterogeneous',
              'Recurrent respiratory infections', 'Fever', 'Hepatomegaly']

sns.barplot(named_keys, filtered_counts)
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('Figures/cftr_shortest_paths_node_counts.png')
plt.show()

counts_df = pd.DataFrame({'node': keys, 'counts': [node_counts[k] for k in keys]})
counts_df['degree'] = [G.degree(k) for k in keys]
counts_df['count:degree ratio'] = counts_df.counts / counts_df.degree
counts_df['name'] = [hpo_g.nodes[k]['name'] if k in hpo_g.nodes() else k for k in keys]
plotting_counts_df = counts_df.sort_values('count:degree ratio', ascending=False).iloc[2:20, :]

sns.barplot(plotting_counts_df['name'], plotting_counts_df['count:degree ratio'])
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('Figures/cftr_shortest_paths_node_counts_ratio.png')
plt.show()

nx.write_edgelist(G, 'Edgelists/string_hpo_CFTR_direct_HPO_connections_removed.txt')
df = pd.read_csv('Edgelists/string_hpo_CFTR_direct_HPO_connections_removed.txt', sep=' ', header=None)
df.columns = ['source', 'target', 'data']
df['predicate'] = '.'
df = df[['source', 'predicate', 'target']]
df.to_csv('Edgelists/string_hpo_CFTR_direct_HPO_connections_removed_triple_list.txt', header=False, index=False,
          sep='\t')

# get depth down HPO
depth = []
for hpo in terms_to_remove:
    depth.append(len(nx.shortest_path(HPO, hpo, 'HP:0000001')))
    nx.degree(g, hpo)


# TODO create a function that takes network, a set of seeds and a set of targets and produces a set of shortest paths a plot of the count / degree ratio of top 20 nodes


def get_shortest_paths(G, seeds, targets):
    paths = []
    # length of paths getting from those terms to CFTR now
    for hpo in seeds:
        for target in targets:
            try:
                temp_paths = list(nx.all_shortest_paths(G, hpo, target))
            except:
                print('No path from ' + hpo + ' to ' + target)
                continue
            if len(list(temp_paths)[0]) > 3:
                print(hpo)
            # print(temp_paths)
            for path in temp_paths:
                # print(len(path))
                # show the visualization
                paths.append(path)

    nodes = []
    node_counts = {}
    for path in paths:
        for p in path:
            nodes.append(p)
            if p in node_counts:
                node_counts[p] += 1
            else:
                node_counts[p] = 1
    return paths, node_counts


def plot_count_degree_node_ratio(counts, fig_name, start_pos=0, end_pos=20):
    """
    :param counts: dictionary keys are node names and values are interger values
        (the node_count value from get_shortest_paths())
    :param fig_name: path/name of output plot
    :param start_pos: integer, start position for ranked node to include in plot
    :param end_pos: integer, end position for ranked node to include in plot
    :return: N/A
    """
    keys = counts.keys()
    counts_df = pd.DataFrame({'node': keys, 'counts': [counts[k] for k in keys]})
    counts_df['degree'] = [G.degree(k) for k in keys]
    counts_df['count:degree ratio'] = counts_df.counts / counts_df.degree
    counts_df['name'] = [hpo_g.nodes[k]['name'] if k in hpo_g.nodes() else k for k in keys]
    plotting_counts_df = counts_df.sort_values('count:degree ratio', ascending=False).iloc[start_pos:end_pos, :]

    sns.barplot(plotting_counts_df['name'], plotting_counts_df['count:degree ratio'])
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(fig_name)
    plt.show()


def get_count_degree_node_ratio_df(counts, start_pos=0, end_pos=20):
    """
    :param counts: dictionary keys are node names and values are interger values
        (the node_count value from get_shortest_paths())
    :return: Pandas dataframe with count:degree ratio
    """
    keys = counts.keys()
    counts_df = pd.DataFrame({'node': keys, 'counts': [counts[k] for k in keys]})
    counts_df['degree'] = [G.degree(k) for k in keys]
    counts_df['count:degree ratio'] = counts_df.counts / counts_df.degree
    counts_df['name'] = [hpo_g.nodes[k]['name'] if k in hpo_g.nodes() else k for k in keys]
    plotting_counts_df = counts_df.sort_values('count:degree ratio', ascending=False).iloc[start_pos:end_pos, :]

    return plotting_counts_df


# CFTR - CFTR
# ps, ncs = get_shortest_paths(G, terms_to_remove, ['CFTR'])
# plot_count_degree_node_ratio(ncs, 'Figures/string_hpo_CFTR-CFTR.png', 2, 20)

# TODO get seeds and targets for HTT and TOF
gene_sets = '2_disease.txt'
disease_gene_sets = {}
for line in open(gene_sets, 'r'):
    row = line.strip().split('\t')
    # the first item is the disease common name, everything else is genes
    disease_gene_sets[row[0]] = row[1:]

aa_seed, aa_target = partition(disease_gene_sets['APLASTIC ANEMIA'], 2)
tof_seed, tof_target = partition(disease_gene_sets['TETRALOGY OF FALLOT'], 2)

htt_hpos = ["HP:0000738", "HP:0000752", "HP:0003324", "HP:0002119", "HP:0011448", "HP:0006855", "HP:0031473",
            "HP:0001336", "HP:0002529", "HP:0002376", "HP:0002340", "HP:0000716", "HP:0003487", "HP:0045082",
            "HP:0000006", "HP:0002073", "HP:0002375", "HP:0001263", "HP:0000007", "HP:0000496", "HP:0002317",
            "HP:0002591", "HP:0002066", "HP:0002312", "HP:0000718", "HP:0002355", "HP:0002300", "HP:0001262",
            "HP:0008936", "HP:0200136", "HP:0000708", "HP:0002059", "HP:0002808", "HP:0003763", "HP:0000741",
            "HP:0001824", "HP:0000734", "HP:0010864", "HP:0001773", "HP:0001288", "HP:0002540", "HP:0002650",
            "HP:0007256", "HP:0002169", "HP:0200147", "HP:0002141", "HP:0000746", "HP:0000737", "HP:0031845",
            "HP:0002354"]

name_index = ['AA', 'TOF', 'HTT', 'CFTR']

if os.path.exists('all_paths.pickle') and os.path.exists('all_node_counts.pickle'):
    all_paths = pickle.load(open('all_paths.pickle', 'rb'))
    all_node_counts = pickle.load(open('all_node_counts.pickle', 'rb'))
else:
    all_paths = {}
    all_node_counts = {}
    for i, seed in enumerate([aa_seed, tof_seed, htt_hpos, terms_to_remove]):
        for j, target in enumerate([aa_target, tof_target, ['HTT'], ['CFTR']]):
            ps, ncs = get_shortest_paths(G, seed, target)
            s_name = name_index[i]
            t_name = name_index[j]
            all_paths[s_name + '-' + t_name] = ps
            all_node_counts[s_name + '-' + t_name] = ncs
            plot_count_degree_node_ratio(ncs, 'Figures/string_hpo_' + s_name + '-' + t_name + '.png', 0, 20)
        pickle.dump(all_node_counts, open('all_node_counts.pickle', 'wb'))
        pickle.dump(all_paths, open('all_paths.pickle', 'wb'))

shortest_paths = np.zeros((4, 4))
keys = list(all_paths.keys())
keys.sort()
for key in keys:
    sets = key.split('-')
    index_1 = name_index.index(sets[0])
    index_2 = name_index.index(sets[1])
    # get average path length
    avg = sum([len(p) for p in all_paths[key]]) / len(all_paths[key])
    shortest_paths[index_1, index_2] = avg
    shortest_paths[index_2, index_1] = avg

shortest_paths_df = pd.DataFrame(shortest_paths)
shortest_paths_df.columns = name_index
shortest_paths_df.index = name_index

plt.figure(figsize=(6, 6))
g = sns.heatmap(shortest_paths_df, annot=True, cbar=False)
plt.xlabel('Seed')
plt.ylabel('Target')
plt.title('Avg Shortest Path Distance')
plt.savefig('Figures/4x4_shuffles.png')
plt.show()

# get the overlap size of the top 20 in CFTR and it's shuffles
cftr_df = get_count_degree_node_ratio_df(all_node_counts['CFTR-CFTR'])
for thing in name_index:
    print(thing)
    df = get_count_degree_node_ratio_df(all_node_counts['CFTR-' + thing], 0, 10000)
    print(len(set(cftr_df.name).intersection(set(df.iloc[0:20, -1]))))
