import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np


sdf = pd.read_csv('all_string_results.csv')
pdf = pd.read_csv('all_pkl_results.csv')

sdf['network'] = 'String'
pdf['network'] = 'PheKnowLater'

df = pd.concat([sdf,pdf])

# subset the columns
percent_df = df[['500 %','100 %', '50 %', '25 %', '10 %', 'disease', 'network']]

# make data tidy
tidy_percent_df = pd.melt(percent_df,['disease','network'], var_name='Top X Genes',value_name="value")
tidy_percent_df['Top X Genes'] = [int(x.replace(' %','') )for x in tidy_percent_df['Top X Genes']]
splot = sns.boxplot(x="network", y="value", hue="Top X Genes", data=tidy_percent_df, palette="Set1")
plt.ylabel('Portion of genes recovered')
plt.savefig('Figures/boxplot.png')
plt.show()

percent_df_2 = percent_df[['500 %','100 %', '50 %', '25 %', '10 %','network']]

r = [1,2]
barWidth = 0.85

mean_series = percent_df_2.groupby('network').agg(pd.Series.mean)
p500 = np.array(list(mean_series['500 %']))
p100 = np.array(list(mean_series['100 %']))
p50 = np.array(list(mean_series['50 %']))
p25 = np.array(list(mean_series['25 %']))
p10 = np.array(list(mean_series['10 %']))

plt.bar(r, p500, color='#df7f20', edgecolor='white', width=barWidth, label='500')
plt.bar(r, p100, bottom=p500, color='#905998', edgecolor='white', width=barWidth, label='100')
plt.bar(r, p50, color='#59a257',bottom=p500+p100, edgecolor='white', width=barWidth, label='50')
plt.bar(r, p25, color='#477ca7',bottom=p500+p100+p50, edgecolor='white', width=barWidth, label='25')
plt.bar(r, p10, color='#cb3336',bottom=p500+p100+p50+p25, edgecolor='white', width=barWidth, label='10')
plt.xticks(r, ['PheKnowLater', 'String'])
plt.xlabel("Network")
plt.ylabel("Average portion of genes rediscovered")
plt.legend(loc='upper left', bbox_to_anchor=(1, 1), ncol=1)
plt.tight_layout()
plt.savefig('Figures/stackedbarplot.png')
plt.show()
