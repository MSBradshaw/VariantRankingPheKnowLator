import sys
import networkx as nx
import random
import re
import pandas as pd
import numpy as np

sys.path.append('pyrwr/')
from pyrwr.rwr import RWR


# https://github.com/jinhongjung/pyrwr

# randomly split a list into n list
def partition(list_in, n):
    random.shuffle(list_in)
    return [list_in[i::n] for i in range(n)]

# my version of random walk... nothing special
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


def create_numbered_edgelist(edge_list, output_file_name):
    nodes = set()
    for line in open(edge_list, 'r'):
        row = line.strip().split('\t')
        nodes.add(row[0])
        nodes.add(row[1])
    nodes = list(nodes)
    nodes.sort()
    # map the node names to a number
    node_mapping = {x: i for i, x in enumerate(nodes)}
    # write it to a file
    with open(output_file_name, 'w') as out:
        for line in open(edge_list, 'r'):
            row = line.strip().split('\t')
            out.write(str(node_mapping[row[0]]))
            out.write('\t')
            out.write(str(node_mapping[row[1]]))
            out.write('\n')
    return node_mapping


def run_rwr(edge_list, start_nodes, target_nodes, node_mapping):
    # node_mapping = create_numbered_edgelist(edge_list)
    num_to_node_mapping = {node_mapping[x]: x for x in node_mapping.keys()}
    start_nums = [node_mapping[x] for x in start_nodes if x in node_mapping]
    target_nums = [node_mapping[x] for x in target_nodes if x in node_mapping]
    # np array to track the sum of residuals from iteration of the loop
    residuals = np.zeros(len(node_mapping))
    rwr = RWR()
    rwr.read_graph(edge_list, 'undirected')
    for node_num in start_nums:
        r = rwr.compute(node_num, c=.15, max_iters=100)
        residuals += r
    r_df = pd.DataFrame({'node': list(range(len(residuals))), 'residuals': residuals})
    r_df = r_df.sort_values('residuals', ascending=False)

    # remove non-genes from rankings
    r_df['node_name'] = [num_to_node_mapping[x] for x in r_df['node']]
    return r_df


def score_network(edgelist, gene_sets, results_name, is_pheknowlater=False):
    disease_gene_sets = {}
    for line in open(gene_sets, 'r'):
        row = line.strip().split('\t')
        # the first item is the disease common name, everything else is genes
        disease_gene_sets[row[0]] = row[1:]

    node_mapping = create_numbered_edgelist(edgelist,edgelist + 'temp_numbered_rwr_edgelist.tsv')
    scores = {'500 count': [], '500 %': [], '100 count': [], '100 %': [], '50 count': [], '50 %': [], '25 count': [],
              '25 %': [], '10 count': [], '10 %': [], 'disease': []}
    num_iterations = 10
    results_df = None
    for disease in disease_gene_sets.keys():
        print(disease)
        for iteration in range(num_iterations):
            print('\t', str(iteration))
            disease_genes = disease_gene_sets[disease]
            nodes = partition(disease_genes, 2)
            start, targets = nodes[0], nodes[1]
            r_df = None
            r_df = run_rwr(edgelist + 'temp_numbered_rwr_edgelist.tsv', start, targets, node_mapping)

            # Pheknowlater needs the non-gene nodes removed
            if is_pheknowlater:
                r_df['is_gene'] = [re.match('.*http.*', x) is None for x in r_df['node_name']]
                r_df = r_df[r_df['is_gene'] == True]

            # the top X we want scores for
            top_xs = [500, 100, 50, 25, 10, 0]
            for i in range(len(top_xs[:-1])):
                top_overlap = len(list(set(r_df.iloc[top_xs[i + 1]:top_xs[i], 2]) & set(targets)))
                scores[str(top_xs[i]) + ' count'].append((top_overlap))
                scores[str(top_xs[i]) + ' %'].append(top_overlap / len(targets))
            scores['disease'].append(disease)
        single_df = pd.DataFrame(scores)
        # get the mean of each of the iterations for this disease
        means = single_df.mean(axis=0)
        scores_mean_df = pd.DataFrame({col: [means[col]] for col in list(single_df.columns)[:-1]})
        scores_mean_df['disease'] = [disease]
        # add the means to the results dataframe
        if results_df is None:
            results_df = scores_mean_df
        else:
            results_df = pd.concat([results_df, scores_mean_df])
    results_df.to_csv(results_name)


# 'pheknowlater_edgelist_gene_symbols_no_snps_clean.tsv'
# 'DisGeNET_genesets.txt'
# 'pheknowlaterresults.csv'

score_network(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]=='True')

# score_network('pheknowlater_edgelist_gene_symbols_no_snps_clean.tsv', 'DisGeNET_genesets.txt', 'pheknowlater_results.csv', True)
# score_network('string_edge_list_common_names.tsv', 'DisGeNET_genesets.txt', 'string_results.csv')

# string_edge_list_common_names.tsv DisGeNET_genesets.txt string_results.csv False
