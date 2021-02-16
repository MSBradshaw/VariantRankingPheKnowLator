import random

gene_sets = 'DisGeNET_genesets.txt'
disease_gene_sets = {}
for line in open(gene_sets, 'r'):
    row = line.strip().split('\t')
    # the first item is the disease common name, everything else is genes
    disease_gene_sets[row[0]] = row[1:]

for key in disease_gene_sets.keys():
    if len(disease_gene_sets[key]) <= 80:
        print(key + '\t' + str(len(disease_gene_sets[key])))

# smallest group that is not cancer or autoimmune: sIMB 'Inclusion Body Myopathy, Sporadic'

# TODO make a shuffled set to test against
all_genes = set()
for key in disease_gene_sets.keys():
    for item in disease_gene_sets[key]:
        if item not in disease_gene_sets['Inclusion Body Myopathy, Sporadic']:
            all_genes.add(item)

shuffled_genes = random.sample(all_genes, len(disease_gene_sets['Inclusion Body Myopathy, Sporadic']))
sIMB_genes = list(disease_gene_sets['Inclusion Body Myopathy, Sporadic'])
with open('sIMB_and_random_shuffle_disease_set.tsv','w') as out:
    out.write('Inclusion Body Myopathy, Sporadic\t')
    for i, item in enumerate(sIMB_genes):
        out.write(item)
        if i != len(sIMB_genes)-1:
            out.write('\t')
    out.write('\n')
    out.write('Randomly Chosen Null Set\t')
    for i, item in enumerate(shuffled_genes):
        out.write(item)
        if i != len(shuffled_genes) - 1:
            out.write('\t')

