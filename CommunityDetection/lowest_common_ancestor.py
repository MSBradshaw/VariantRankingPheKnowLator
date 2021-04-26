import networkx as nx
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import obonet


def get_coms(file):
    coms = {}
    for line in open(file, 'r'):
        row = line.strip().split('\t')
        for com in row[1:]:
            if com in coms:
                coms[com].append(row[0])
            else:
                coms[com] = [row[0]]
    return coms

def to_coms_file(coms,outfile_name):
    with open(outfile_name,'w') as file:
        for x in coms:
            file.write(x + '\t')
            file.write('\t'.join(coms[x]))
            file.write('\n')

def to_single_com_file(coms,outfile_base_name):
    for x in coms:
        with open(outfile_base_name + str(x) + '.txt','w') as file:
            for y in coms[x]:
                file.write(y + '\n')

m_coms = get_coms('CommunityDetection/markov.txt')
l_coms = get_coms('CommunityDetection/louvain.txt')
w_coms = get_coms('CommunityDetection/walktrap.txt')

url = 'https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.obo'
hpo_g = obonet.read_obo(url)

hpos = ['HP:0002865', 'HP:0002890']


def lowest_common_ancestor(hpos, G):
    winning_count = len(hpos)
    ancestors = {}
    currents = hpos.copy()
    while True:
        anc = sum(ancestors[k] if ancestors[k] >= winning_count else 0 for k in ancestors.keys())
        if anc != 0:
            # get the name of the lca
            temp = [k for k in ancestors.keys() if ancestors[k] >= winning_count]
            return temp[0]

        for i, h in enumerate(hpos):
            try:
                currents[i] = G.nodes[currents[i]]['is_a'][0]
            except KeyError:
                winning_count -= 1
            if currents[i] not in ancestors:
                ancestors[currents[i]] = 1
            else:
                ancestors[currents[i]] += 1


def get_average_lca(coms, G,outname,title):
    lcas = []
    for k in coms.keys():
        com = coms[k]
        if sum('HP:' in x for x in com) >= 2:
            parent = lowest_common_ancestor(com, G)
            for x in com:
                if 'HP:' not in x:
                    continue
                try:
                    # print(len(nx.shortest_path(G, x, parent)))
                    lcas.append(len(nx.shortest_path(G, x, parent)))
                except nx.exception.NodeNotFound:
                    pass
    plt.hist(lcas)
    plt.xlabel('Average Distance Up HPO')
    plt.ylabel('Count')
    plt.title(title)
    plt.savefig(outname)
    plt.show()
    return sum(lcas) / len(lcas)

gg = nx.DiGraph.to_undirected(hpo_g)
print('Markov', get_average_lca(m_coms, gg,'Figures/m_lva_avg.png','Markov'))
print('Louvain', get_average_lca(l_coms, gg,'Figures/l_lva_avg.png','Louvain'))
print('Walktrap', get_average_lca(w_coms, gg,'Figures/w_lva_avg.png','Walktrap'))
lowest_common_ancestor(['HP:0002865', 'HP:0002890'], hpo_g)
