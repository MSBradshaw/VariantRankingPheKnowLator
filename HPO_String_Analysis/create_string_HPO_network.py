import obonet
import pandas as pd
import networkx as nx


url = 'https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.obo'
hpo_g = obonet.read_obo(url)

hpo_to_gene_df = pd.read_csv('HPO_String_Analysis/genes_to_phenotype.txt',sep='\t',comment='#')

hpo_to_gene_df.columns = ['entrez-gene-id','entrez-gene-symbol','HPO-Term-Name','HPO-Term-ID','Frequency-Raw','Frequency-HPO','Additional Info from G-D source','G-D source','disease-ID for link']

for i in range(hpo_to_gene_df.shape[0]):
    hpo_g.add_edge(hpo_to_gene_df.iloc[i,1],hpo_to_gene_df.iloc[i,2])

# add string to the graph
for line in open('Edgelists/string_edge_list_common_names.tsv','r'):
    row = line.strip().split('\t')
    hpo_g.add_edge(row[0],row[1])

nx.write_edgelist(hpo_g,'HPO_String_Analysis/HPO_String_edgelist.tsv',delimiter='\t')
