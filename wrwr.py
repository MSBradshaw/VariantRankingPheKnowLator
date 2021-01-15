import sys
import networkx as nx
import random
import re
import pandas as pd
import numpy as np
import pickle
from numpy.random import choice
import time
import argparse


# randomly split a list into n list
def partition(list_in, n):
    random.shuffle(list_in)
    return [list_in[i::n] for i in range(n)]

def append_path_to_file(path,file):
    with open(file,'a') as out:
        out.write(','.join(path))
        out.write('\n')

# my version of random walk... nothing special
def wrwr(G, seeds, reset_out_of, threshold):
    walk_res = {n: 1 for n in G.nodes}
    previous_weights = walk_res.copy()
    previous_sum_of_weights = len(G.nodes)
    iterations = 0
    current = None
    path = []
    while True:
        if random.randint(0, reset_out_of) == 0 or current is None:
            append_path_to_file(path,str(threshold) + '_wrwr_paths.csv')
            current = random.choice(seeds)
            path = [current]
            # print('Restarting')
            continue
        # choose a random seed node
        try:
            neighs = list(nx.all_neighbors(G, current))
        except nx.exception.NetworkXError:
            current = None
            continue
        weights = [G.edges[current, n]['weight'] for n in neighs]
        sw = sum(v for v in weights)
        weights = [w / sw if sw != 0 else 0 for w in weights]
        # if there are no neighbors, you have reached a leaf in the graph, end the walk.
        if len(neighs) == 0:
            # print('No neighbors')
            current = None
            continue
        current = choice(neighs, 1, weights)[0]
        path.append(current)
        if current not in walk_res:
            walk_res[current] = 1
        else:
            walk_res[current] += 1
        if iterations % 10 == 0:
            sum_of_weights = sum(walk_res[n] for n in walk_res.keys())
            sum_of_diff = sum(
                previous_weights[n] / previous_sum_of_weights - walk_res[n] / sum_of_weights for n in walk_res.keys())
            previous_sum_of_weights = sum_of_weights
            previous_weights = walk_res
            print(sum_of_diff)
            if sum_of_diff < threshold and iterations > 100:
                break
        iterations += 1
    return walk_res


def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--edge_list",
                        dest="edge_list",
                        required=True,
                        default=None,
                        help="path to edge list")

    parser.add_argument("--disease_sets",
                        dest="disease_sets",
                        required=True,
                        default=None,
                        help="path to disease set file")

    parser.add_argument("--edge_weights",
                        dest="edge_weights",
                        required=False,
                        default=None,
                        help="path to edgelist")

    parser.add_argument("--node_mapping",
                        dest="node_mapping",
                        required=False,
                        default=None,
                        help="path to node and egde list name mapping")

    parser.add_argument("--threshold",
                        dest="threshold",
                        required=True,
                        default=None,
                        help="threshold for when the walks will be considered as having converged")


    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    if args.node_mapping:
        url_to_common_name = {line.strip().split()[0]: line.strip().split()[1] for line in open(args.node_mapping)}
    # read in a network with edge weights for certain types 'Edgelists/pkl_triples_with_symbols.tsv'
    G = nx.Graph()
    el = args.edge_list
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
    # read in a set of weights here and re assign all edge weights: 'proportional_change_weights.tsv'
    if args.edge_weights:
        weight_mapping = {line.strip().split('\t')[0]:line.strip().split('\t')[1] for line in open(args.edge_weights,'r')}
        for edge in G.edges:
            if G.edges[edge]['type'] in weight_mapping:
                G.edges[edge]['weight'] = weight_mapping[G.edges[edge]['type']]

    gene_sets = args.disease_sets
    disease_gene_sets = {}
    for line in open(gene_sets, 'r'):
        row = line.strip().split('\t')
        # the first item is the disease common name, everything else is genes
        disease_gene_sets[row[0]] = row[1:]

    ranked_gene_names = {'gene': [], 'rank': [], 'disease': [], 'is_target': []}
    scores = {'500 count': [], '500 %': [], '100 count': [], '100 %': [], '50 count': [], '50 %': [], '25 count': [],
              '25 %': [], '10 count': [], '10 %': [], 'disease': []}
    disease = list(disease_gene_sets.keys())[0]
    disease_genes = disease_gene_sets[list(disease_gene_sets.keys())[0]]
    nodes = partition(disease_genes, 2)
    start, targets = nodes[0], nodes[1]
    good_starts = [s for s in start if s in G.nodes]
    start_time = time.time()
    res = wrwr(G, good_starts, 10, float(args.threshold))
    end = time.time()
    print(end - start_time)
    r_df = pd.DataFrame({'node': list(res.keys()), 'residuals': list(res.values())})
    r_df = r_df.sort_values('residuals', ascending=False)

    ranked_gene_names['gene'] += list(r_df['node'])
    ranked_gene_names['rank'] += list(range(r_df.shape[0]))
    ranked_gene_names['disease'] += [disease] * r_df.shape[0]
    ranked_gene_names['is_target'] += [x in targets for x in r_df['node']]

    # get the number of genes in the disease gene set
    num_targets_in_network = sum(t in list(r_df['node']) for t in targets)

    # the top X we want scores for
    top_xs = [500, 100, 50, 25, 10, 0]
    for i in range(len(top_xs[:-1])):
        print(i)
        top_overlap = len(list(set(r_df.iloc[top_xs[i + 1]:top_xs[i], 0]) & set(targets)))
        scores[str(top_xs[i]) + ' count'].append((top_overlap))
        scores[str(top_xs[i]) + ' %'].append(top_overlap / num_targets_in_network)
    scores['disease'].append(disease)
    scores['time'] = [end - start_time]
    scores['threshold'] = [float(args.threshold)]
    pd.DataFrame(scores).to_csv('weighted_converging_test'+args.threshold+'.csv')
    # get a set of weights
    # run WRWR to convergence for one set
    # report residuals
    # score
    # 0.000015

    quit()
    # CONVERT Predicates to gene symbols
    id_2_symbol_dict = {line.strip().split('\t')[0]: line.strip().split('\t')[1] for line in
                        open('Data/gene_id2symbol.txt', 'r')}
    with open('Edgelists/pkl_triples_with_symbols.tsv', 'w') as out:
        for line in open(
                '/Users/michael/PycharmProjects/VariantRankingPheKnowLator/Data/PheKnowLator_Subclass_InverseRelations_NotClosed_NoOWL_Triples_Identifiers_11_24_2020.txt',
                'r'):
            if 'www.ncbi.nlm.nih.gov/snp/' in line:
                continue
            row = line.strip().split()
            row[0] = row[0].replace('https://www.ncbi.nlm.nih.gov/gene/', '')
            row[2] = row[2].replace('https://www.ncbi.nlm.nih.gov/gene/', '')
            if row[0] in id_2_symbol_dict:
                row[0] = id_2_symbol_dict[row[0]]
            if row[2] in id_2_symbol_dict:
                row[2] = id_2_symbol_dict[row[2]]
            out.write(row[0] + '\t' + row[1] + '\t' + row[2] + '\n')
