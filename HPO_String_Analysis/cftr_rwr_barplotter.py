import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

ndf = pd.read_csv('HPO_String_Analysis/normal_wrwr_0.000010.csv')
adf = pd.read_csv('HPO_String_Analysis/anti_degree_0.000011.csv')

ndf['network'] = 'Normal'
adf['network'] = 'Anti Degree'

df = pd.concat([ndf, adf])

# subset the columns
percent_df = df[['500 %', '100 %', '50 %', '25 %', '10 %', 'disease', 'network']]

# make data tidy
tidy_percent_df = pd.melt(percent_df, ['disease', 'network'], var_name='Top X Genes', value_name="value")
tidy_percent_df['Top X Genes'] = [int(x.replace(' %', '')) for x in tidy_percent_df['Top X Genes']]
splot = sns.barplot(x="network", y="value", hue="Top X Genes", data=tidy_percent_df, palette="Set1")
plt.ylabel('Portion of genes recovered')
plt.savefig('Figures/cftr_normal_vs_anti_degree_rwr_barplot.png')
plt.show()
