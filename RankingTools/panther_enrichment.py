import requests
import json
import pandas as pd
import argparse


def get_args():
    parser = argparse.ArgumentParser()


    parser.add_argument('--coms',
                        dest='coms',
                        help='tsv file where each row in a community, first area is com name followed by a tab ' +
                             'seporated list of community members (not all rows will be the same length)')

    parser.add_argument('--pval',
                        dest='pval',
                        help='threshold for significance')

    parser.add_argument('--prefix',
                        dest='prefix',
                        help='prefix to be used on each output file')

    args = parser.parse_args()
    return args


def go_enrichment(genes):
    """
    :param genes: list of genes
    :return: pandas DataFrame of the panther API results, not all returned rows are significant
    """
    g_string = ','.join(genes)
    url = 'http://pantherdb.org/services/oai/pantherdb/enrich/overrep?geneInputList=' + g_string +\
               '&organism=9606&annotDataSet=GO%3A0008150&enrichmentTestType=FISHER&correction=FDR'

    resp = requests.get(url)
    resp_obj = json.loads(resp.content)

    results = {'number_in_list': [],
               'fold_enrichment': [],
               'fdr': [],
               'expected': [],
               'number_in_reference': [],
               'pValue': [],
               'id': [],
               'label': [],
               'plus_minus': []}
    try:
        for i in range(len(resp_obj['results']['result'])):
            for key in resp_obj['results']['result'][i].keys():
                if key != 'term':
                    results[key].append(resp_obj['results']['result'][i][key])
                else:
                    try:
                        results['id'].append(resp_obj['results']['result'][i]['term']['id'])
                    except KeyError:
                        results['id'].append('None')
                    try:
                        results['label'].append(resp_obj['results']['result'][i]['term']['label'])
                    except KeyError:
                        results['label'].append('None')

    except KeyError:
        # there are no results, return empty dataframe
        return pd.DataFrame()
    return pd.DataFrame(results)


if __name__ == "__main__":
    args = get_args()
    for line in open(args.coms,'r'):
        row = line.strip().split('\t')
        genes = [x for x in row[1:] if 'HP:' not in x]
        if len(genes) > 5000:
            print('Too many members')
            continue
        df = go_enrichment(genes)
        df=df[df['pValue'] < float(args.pval)]
        df.to_csv(args.prefix + row[0] + '.tsv', sep='\t', index=False)

quit()
"""
# example function use
example_genes = ['SPAG1','TTC12','OFD1','DPCD','FOXJ1','WDR63','POC1B','FGFR1OP','DNAJB13']
example_df = go_enrichment(example_genes)

# example command line use
python panther_enrichment.py --coms del_coms.tsv --prefix GO/Results/panther_test --pval 0.0005
"""
