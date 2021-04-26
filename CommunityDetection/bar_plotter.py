import networkx as nx
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def get_coms(file):
    coms = {}
    for line in open(file, 'r'):
        row = line.strip().split('\t')
        for com in row[1:]:
            if com in coms:
                coms[com].append(row[0])
            else:
                coms[com] = [row[0]]
    return coms


m_coms = get_coms('CommunityDetection/markov.txt')
l_coms = get_coms('CommunityDetection/louvain.txt')
w_coms = get_coms('CommunityDetection/walktrap.txt')

m_sizes = [len(m_coms[k]) for k in m_coms.keys()]
l_sizes = [len(l_coms[k]) for k in l_coms.keys()]
w_sizes = [len(w_coms[k]) for k in w_coms.keys()]

print(len(m_coms['0']))
print(len(l_coms['0']))
print(len(w_coms['0']))

plt.scatter(list(range(len(m_sizes))), m_sizes)
plt.xlabel('community id')
plt.ylabel('# of members')
plt.yscale('log')
plt.title('Number of communities=' + str(len(m_sizes)))
plt.savefig('CommunityDetection/markov_coms.png')
plt.show()

plt.scatter(list(range(len(l_sizes))), l_sizes)
plt.yscale('log')
plt.xlabel('community id')
plt.ylabel('# of members')
plt.title('Number of communities=' + str(len(l_sizes)))
plt.savefig('CommunityDetection/louvain_coms.png')
plt.show()

plt.scatter(list(range(len(w_sizes))), w_sizes)
plt.yscale('log')
plt.xlabel('community id')
plt.ylabel('# of members')
plt.title('Number of communities=' + str(len(w_sizes)))
plt.savefig('CommunityDetection/walktrap_coms.png')
plt.show()

# get the ratio of HPO to protein
m_ratios = [sum('HP:' in x for x in m_coms[com]) / len(m_coms[com]) for com in m_coms]
l_ratios = [sum('HP:' in x for x in l_coms[com]) / len(l_coms[com]) for com in l_coms]
w_ratios = [sum('HP:' in x for x in w_coms[com]) / len(w_coms[com]) for com in w_coms]
box_df = pd.DataFrame({'group':['Markov'] * len(m_ratios) +
                               ['Louvain'] * len(l_ratios) +
                               ['Walktrap'] * len(w_ratios),
                       'HPO/Gene member count ratio': m_ratios + l_ratios + w_ratios})
sns.boxplot('group','HPO/Gene member count ratio',data=box_df)
plt.show()

plt.hist(m_ratios)
plt.xlabel('HPO/Gene member count ratio')
plt.ylabel('Count')
plt.title('Markov')
plt.savefig('CommunityDetection/m_ratio_hist.png')
plt.yscale('log')
plt.savefig('CommunityDetection/m_ratio_hist_log.png')
plt.show()

plt.hist(l_ratios)
plt.xlabel('HPO/Gene member count ratio')
plt.ylabel('Count')
plt.title('Louvain')
plt.savefig('CommunityDetection/l_ratio_hist.png')
plt.yscale('log')
plt.savefig('CommunityDetection/l_ratio_hist_log.png')
plt.show()

plt.hist(w_ratios)
plt.xlabel('HPO/Gene member count ratio')
plt.ylabel('Count')
plt.title('Walktrap')
plt.yscale('log')
plt.savefig('CommunityDetection/w_ratio_hist.png')
plt.yscale('log')
plt.savefig('CommunityDetection/w_ratio_hist_log.png')
plt.show()



def variance(data):
    # Number of observations
    n = len(data)
    # Mean of the data
    mean = sum(data) / n
    # Square deviations
    deviations = [(x - mean) ** 2 for x in data]
    # Variance
    v = sum(deviations) / n
    return v

variance(m_ratios)