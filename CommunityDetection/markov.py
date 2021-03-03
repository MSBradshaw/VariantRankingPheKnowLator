from cdlib import algorithms
import networkx as nx
import pickle

# G = nx.read_edgelist('HPO_String_Analysis/HPO_String_edgelist.tsv', delimiter='\t')
G = nx.karate_club_graph()

coms = algorithms.markov_clustering(G)
pickle.dump(coms, open('CommunityDetection/louvain_coms.pickle', 'wb'))

cd = coms.to_node_community_map()

with open('markov.txt','w') as outfile:
    for key in cd.keys():
        outfile.write(key + '\t' + '\t'.join([str(x) for x in cd[key]]) + '\n')
