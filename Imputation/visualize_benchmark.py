import os 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# font size
plt.rcParams.update({'font.size': 28})

samples = ['JS34', 'JS35', 'JS36', 'JS40']
n_neighbors_list = [5,10,20,25,40,50,100,200,400]

if not os.path.exists('figures/imputation'):
    os.makedirs('figures/imputation')

for sample in samples:
    fig, ax = plt.subplots(figsize = (12, 8))
    for n_neighbors in n_neighbors_list:
        accu_pcc = []
        df_pcc = pd.read_csv(f'/home/jinmr2/sample_integration/three_samples/imputation/v5_seurat_all_genes/output/pearsonrGenes_{sample}_{n_neighbors}.csv', header = 0, index_col = 0)
        df_pcc = df_pcc.sort_values(by = 'pearsonr', ascending = True)
        pcc_scores = df_pcc['pearsonr'].values
        pcc_scores = pcc_scores[~pd.isna(pcc_scores)]
        for i in range(len(pcc_scores)):
            if i == 0:
                accu_pcc.append(pcc_scores[i])
            else:
                accu_pcc.append(accu_pcc[-1] + pcc_scores[i])
        ax.plot(accu_pcc, label = f'n_neighbors = {n_neighbors}')
    ax.set_xlabel('Number of Genes')
    ax.set_ylabel('Imputation Performance Score')

    ax.legend(markerscale = 20, fontsize = 20)
    fig.savefig(f'figures/imputation/imputation_performance_score_{sample}.png', dpi = 300, bbox_inches = 'tight')
