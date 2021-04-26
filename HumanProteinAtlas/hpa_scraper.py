import urllib.request
import pandas as pd
from bs4 import BeautifulSoup


norm_tissues = pd.read_csv('HumanProteinAtlas/normal_tissue.tsv', sep='\t')
url_parts = list(set(norm_tissues['Gene'] + '-' + norm_tissues['Gene name']))
base_url = 'https://www.proteinatlas.org/'

genes = []
tissue_specificity = []
for part in url_parts:
    try:
        genes.append(part.split('-')[1])
        url = base_url + part
        uf = urllib.request.urlopen(url)
        html = uf.read()
        soup = BeautifulSoup(html, 'lxml')
        # the second table, 4th row, first value, first div's text (how to find the tissue specificity table row)
        ts = soup.find_all('table')[1].find_all('tr')[3].find_all('td')[0].find_all('div')[0].text
        tissue_specificity.append(ts)
    except:
        print(part)

pd.DataFrame({'Gene': genes, 'Tissue Specificity': tissue_specificity}).to_csv('tissue_specificity.csv', sep='\t',
                                                                               index=False)
with open('gene_names.txt','w') as file:
    file.write(','.join(list(set(norm_tissues['Gene name']))))