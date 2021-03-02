from cdlib import algorithms
import networkx as nx
import pickle

G = nx.read_edgelist('HPO_String_Analysis/HPO_String_edgelist.tsv', delimiter='\t')
# G = nx.karate_club_graph()

coms = algorithms.louvain(G, weight='weight', resolution=1., randomize=False)
pickle.dump(coms, open('CommunityDetection/louvain_coms.pickle', 'wb'))
