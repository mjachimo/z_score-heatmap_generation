import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import cm
from adjustText import adjust_text
import statistics as st

# Extracts and z-score normalizes significantly enriched normalized counts

df = pd.read_csv('F-vs-M-all.gene_formated.csv')
significance = df['q-value'].tolist()
gene_name = df['gene_id'].tolist()

counts_df = pd.read_csv('counts_anno.csv')
counts_gene_name = counts_df['id'].tolist()

# extracts column names with count data
header = counts_df.columns.values
counts_headers = []
for x in range(len(header)):
    if 'counts' in header[x]:
        counts_headers.append(header[x])
new_df = counts_df[counts_headers]
print(counts_headers)
important_genes = []
new_counts_headers = [counts_headers[3], counts_headers[4], counts_headers[5], counts_headers[0], counts_headers[1], counts_headers[2]]
final_df = pd.DataFrame(columns=new_counts_headers)

# extracts only significant gene counts and z-score normalizes them
for x in range(len(gene_name)):
    if significance[x] < 0.05:
        for k in range(len(new_df)):
            if gene_name[x] == counts_gene_name[k]:
                important_genes.append(gene_name[x])
                gene_counts_values = []
                normalized_counts_values = []
                for y in range(6):
                    gene_counts_values.append(new_df.iloc[k, y])
                mean = st.mean(gene_counts_values)
                std = st.stdev(gene_counts_values)
                for t in range(6):
                    normalized_counts_values.append((new_df.iloc[k, t] - mean) / std)
        new_row = {counts_headers[3]: normalized_counts_values[3], counts_headers[4]: normalized_counts_values[4], counts_headers[5]: normalized_counts_values[5], counts_headers[0]: normalized_counts_values[0], counts_headers[1]: normalized_counts_values[1],
            counts_headers[2]: normalized_counts_values[2]}
        final_df = final_df.append(new_row, ignore_index=True)

print(len(gene_name))
print(len(important_genes))
df2 = pd.DataFrame(important_genes, columns=['gene_names'])
df2.to_csv('gene_names_2.csv', index=False)
final_df.to_csv('normalized expression_2.csv', index=False)


# plots the normalized expression on a heatmap
new_df = pd.read_csv('normalized expression_2.csv')
gene_df = pd.read_csv('gene_names_2.csv')
gene_names = gene_df['gene_names'].tolist()
print(new_df)
of_interest = pd.read_csv('of_interest.csv')
labels = of_interest['gene_names'].tolist()

fig, ax = plt.subplots(figsize=(10, 10))
sns.heatmap(new_df, ax=ax, annot=False, xticklabels=False, yticklabels=False, cmap=sns.color_palette("RdBu_r", as_cmap=True), cbar_kws={'shrink': 0.3,'aspect': 10, 'location': 'left', 'label': 'Normalized Gene counts', 'ticks':[-2, -1, 0, 1, 2]})
z = []
order = []
for x in range(len(labels)):
    order.append(x)
    tally = 0
    for y in range(len(gene_names)):
        if labels[x] == gene_names[y]:
            z.append(y)
        else:
            tally = tally + 1
    if tally == len(gene_names):
        print(labels[x])

z, labels = (list(t) for t in zip(*sorted(zip(z, labels))))
spacer = 15
cell_height = z.copy()
problem = False
x = 0

while x < len(z):
    y = 0
    while y < len(z):
        if ((z[x] - spacer <= z[y] <= z[x] + spacer) and (x != y)):
            z[y] = z[y] + 15
            x = 0
            print(z)
            break
        y = y + 1
    x = x + 1
ax = [ax.annotate(labels[i], (6, cell_height[i]), xytext=(6.5, z[i]), arrowprops=dict(arrowstyle='-', lw=1), ha= 'left', va='center', size=8) for i, txt in enumerate(z)]

plt.savefig('F_vs_M_heatmap_2.png')
plt.show()
