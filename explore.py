import networkx as nx
import pickle
import random
import rdflib
import re
import pandas as pd

df = pd.read_csv('Data/PheKnowLator_Subclass_RelationsOnly_NotClosed_NoOWL_Triples_Identifiers.txt', sep='\t')
df.pop('predicate')
df.to_csv('edgelist.tsv', sep='\t', index=False, header=False)

# remove snps
with open('edgelist_snp_free.tsv', 'w') as file:
    for line in open('edgelist.tsv', 'r'):
        if not re.match('.*\/snp\/.*', line) and not re.match('.*http://purl.obolibrary.org/obo/BFO_0000001.*',
                                                              line) and not re.match(
            '.*http://purl.obolibrary.org/obo/doid.*', line):
            file.write(line)

# re.match('.*\/snp\/.*', line) and
# not re.match('.*http://purl.obolibrary.org/obo/BFO_0000001.*',line)
# read in network as digraph
G = nx.read_edgelist('edgelist_snp_free.tsv', delimiter='\t', create_using=nx.DiGraph())

# dictionary for storing results
walk_res = {}

# starting position for RWR
K = 10
start_node = 'https://www.ncbi.nlm.nih.gov/gene/91624'
start_nodes = ['https://www.ncbi.nlm.nih.gov/gene/91624', 'http://purl.obolibrary.org/obo/GO_0046632']
for seed_node in start_nodes:
    current = seed_node
    for i in range(10):
        count = 0
        while True:
            if count == K:
                break
            count += 1
            # print(current)
            # do the walk, choose one neighbor at random
            neighs = list(nx.all_neighbors(G, current))
            # if there are no neighbors, you have reached a leaf in the graph, end the walk.
            if len(neighs) == 0:
                print('No neighbors')
                break
            current = random.choice(neighs)
            if current not in walk_res:
                walk_res[current] = 1
            else:
                walk_res[current] += 1
            # the random chance of reset
            if random.randint(0, 10) == 0:
                current = seed_node
                print('Restart!')


# make an RWR function
# run RWR for each input gene
# rank the genes
# plot the subgraph explored in RWR with node weights

# do not visit nodes twice in the same walk

def random_walk_with_reset(G, seeds, num_iterations, reset_out_of, K):
    walk_res = {}
    for seed_node in seeds:
        current = seed_node
        for i in range(num_iterations):
            count = 0
            while True:
                if count == K:
                    break
                count += 1
                # print(current)
                # do the walk, choose one neighbor at random
                neighs = list(nx.all_neighbors(G, current))
                # if there are no neighbors, you have reached a leaf in the graph, end the walk.
                if len(neighs) == 0:
                    print('No neighbors')
                    break
                current = random.choice(neighs)
                if current not in walk_res:
                    walk_res[current] = 1
                else:
                    walk_res[current] += 1
                # the random chance of reset
                if random.randint(0, reset_out_of) == 0:
                    current = seed_node
                    print('Restart!')
    return walk_res


r = random_walk_with_reset(G, start_nodes, 10, 10, 10)

# check out the walker from Ideker's lab
# https://github.com/idekerlab/Network_Evaluation_Tools/blob/master/Network%20Evaluation%20Examples/run_network_evaluation.py

genes = set(re.sub('https://www.ncbi.nlm.nih.gov/gene/', '', n) + '\n' for n in G.nodes if re.match('.*gene.*', n))
with open('gene_urls.txt', 'w') as out:
    for n in genes:
        out.write(n + '\n')

# that file was then converted symbols using the web interface https://syngoportal.org/convert.html

id_symbol_map = {line.strip().split('\t')[0]: line.strip().split('\t')[1] for line in open('gene_id2symbol.txt', 'r')}

f = 0
nf = 0
nf_set = set()
# rename nodes to use gene symbol rather than id
with open('pheknowlater_edgelist_gene_symbols_no_snps.tsv', 'w') as out:
    for line in open('edgelist_snp_free.tsv'):
        line = line.strip()
        edge = line.split('\t')
        outline = ''
        if re.match('.*gene.*', edge[0]):
            if re.sub('https://www.ncbi.nlm.nih.gov/gene/', '', edge[0]) in id_symbol_map:
                outline += id_symbol_map[re.sub('https://www.ncbi.nlm.nih.gov/gene/', '', edge[0])] + '\t'
                f += 1
            else:
                nf += 1
                nf_set.add(edge[0])
                # print('Not found ', edge[0])
        else:
            outline += edge[1] + '\t'

        if re.match('.*gene.*', edge[1]):
            if re.sub('https://www.ncbi.nlm.nih.gov/gene/', '', edge[1]) in id_symbol_map:
                outline += id_symbol_map[re.sub('https://www.ncbi.nlm.nih.gov/gene/', '', edge[1])] + '\n'
                f += 1
            else:
                nf += 1
                nf_set.add(edge[1])
                # print('Not found ', edge[1])
        else:
            outline += edge[1] + '\n'
        out.write(outline)

# TODO find out why there are not tripples in the edge list
with open('pheknowlater_edgelist_gene_symbols_no_snps_clean.tsv', 'w') as out:
    for line in open('pheknowlater_edgelist_gene_symbols_no_snps.tsv'):
        row = line.strip().split('\t')
        if len(row) == 2:
            out.write(line)

# convert nodes from names to numbers
nodes = set()
for line in open('pheknowlater_edgelist_gene_symbols_no_snps_clean.tsv', 'r'):
    row = line.strip().split('\t')
    nodes.add(row[0])
    nodes.add(row[1])

nodes = list(nodes)
nodes.sort()
# map the node names to a number
node_mapping = {x:i for i,x in enumerate(nodes)}
pickle.dump(node_mapping,open('pheknowlater_node_number_mapping.pickle','wb'))
# write to file with node numbers
with open('pheknowlater_node_number_edgelist.tsv', 'w') as out:
    for line in open('pheknowlater_edgelist_gene_symbols_no_snps_clean.tsv', 'r'):
        row = line.strip().split('\t')
        out.write(str(node_mapping[row[0]]))
        out.write('\t')
        out.write(str(node_mapping[row[1]]))
        out.write('\n')