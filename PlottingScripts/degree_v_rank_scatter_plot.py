import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from scipy import stats

string = pd.read_csv('Results/string_ranked_results.csv')
string = string[string['rank'] <= 500]
# remove genes that are not in the top 500

plk = pd.read_csv('Results/pkl_rank_results.csv')
plk = plk[plk['rank'] <= 500]

plk_G = nx.read_edgelist('Edgelists/pheknowlater_edgelist_gene_symbols_no_snps_clean.tsv')
string_G = nx.read_edgelist('Edgelists/string_edge_list_common_names.tsv')

string['degree'] = [ string_G.degree[x] for x in string['gene']]
plk['degree'] = [ plk_G.degree[x] for x in plk['gene']]
# get the degree for each gene in the get work

plt.scatter(string['rank'],string['degree'])
plt.ylabel('Degree')
plt.xlabel('Rank')
plt.title('Gene By Rank and Gene in StringDB')
plt.savefig('Figures/scatter_degree_rank_string.png')
plt.show()

plt.scatter(plk['rank'],plk['degree'])
plt.ylabel('Degree')
plt.xlabel('Rank')
plt.title('Gene By Rank and Gene in PheKnowLater')
plt.savefig('Figures/scatter_degree_rank_plk.png')
plt.show()

stats.pearsonr(plk['rank'],plk['degree'])
stats.pearsonr(string['rank'],string['degree'])

# TODO find disease gene set that worked well and one that did not work
