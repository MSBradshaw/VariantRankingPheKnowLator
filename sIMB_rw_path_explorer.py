import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

path_counts = {}
paths = []
lengths = []
with open('1e-05_wrwr_paths.csv','r') as file:
    for line in file:
        line=line.strip()
        if line in path_counts:
            path_counts[line] += 1
        else:
            path_counts[line] = 1
        paths.append(line.split(','))
        lengths.append(len(paths[-1]))

# path length
plt.boxplot(lengths)
plt.ylabel('path length')
plt.savefig('Figures/sIMB_rw_path_lengths_boxplot2.png')
plt.show()

# most common paths
path_counts_df = pd.DataFrame({'path':list(path_counts.keys()),'count':list(path_counts.values())}).sort_values('count',ascending=False)
path_counts_df.iloc[21:40,:]
