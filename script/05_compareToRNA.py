#!/usr/bin/env python
### Comparing between Genes with ATAC-seq peak and genes expression


import pandas as pd
import scipy.stats as ss
import numpy as np
import random
import subprocess, sys
import time
import argparse
import subprocess
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
from matplotlib.ticker import FuncFormatter




def main():
        parser = argparse.ArgumentParser()
        parser.add_argument('atac_gene',help='ATAC-seq peak containing genes in bed',type=str)
        parser.add_argument('RNA_seq',help='Expression table',type=str)
        parser.add_argument('outVenn',help='Venn diagram output file name',type=str)
        parser.add_argument('outBar',help='Barplot output file name',type=str)
        parser.add_argument('-b','--bins',help='Number of bins; default 10',type=int,default=10)

        args = parser.parse_args()

        A = args.atac_gene
        R = args.RNA_seq
        B = args.bins
        O2 = args.outVenn
        O = args.outBar

#        A = 'hESC_Ctrl_1_peak_gene_list.txt'
#        R = 'hESC_Ctrl_1_rpkm.txt'
#        B = 10
#        O = 'Barplot_ATACseq_RNA.png'
        atac=pd.read_csv(A,sep="\t",header=None)
        rna=pd.read_csv(R,sep="\t",header=None)
        rna = rna.sort_values(1).reindex()
        num_gene = rna[0].count()
#        rna[2] = np.sort([random.randint(0,10) for i in range(rna[0].count())])
        rna[2] = np.sort(list(range(B)) * (int(num_gene / B) + 1))[0:num_gene]
        rna[3] = [1 if item in atac[3].to_list() else 0 for item in rna[0]]
        count = [rna[rna[2] == x][3].sum()  for x in range(B)]
        size = [rna[rna[2] == x][3].count()  for x in range(B)]
        ratio = [rna[rna[2] == x][3].sum()/1.0/rna[rna[2] == x][3].count()  for x in range(B)]
        x = np.array(range(B))

        plt.style.use('classic')
        fig = plt.figure()
        ax=fig.add_subplot(1,1,1)
        loc = [str(y) for y in range(1,B+1)]
        ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _:'{:.0%}'.format(y)))
        ax.bar(x, ratio, color='royalblue',align='center')
        plt.xlabel('Expression groups (low to high)')
        plt.ylabel('Percentage of accessible genes')
        plt.xlim(-1,B)
        plt.ylim(0,1)
        ax.set_xticks(x)
        ax.set_xticklabels(loc)
        ax.set_title('ATAC-seq vs. RNA-seq', fontsize=20)
        plt.savefig(O,dpi=300)
        plt.close(fig)

        atac_count = atac[0].count()
        venns_value = [size[0]-count[0],size[B-1]-count[B-1],0,atac_count - count[0] - count[B-1],count[0],count[B-1],0]
        figV = plt.figure(figsize=(4,4))
        v = venn3(subsets = venns_value, set_labels = ('Lowly expressed genes', 'Highly expressed genes', 'Accessible (open) genes'),set_colors=('skyblue','lightgreen','orange'), alpha = 0.7)
        plt.title("Chromatin accessibility vs. Expression")
        for text in v.set_labels:
            text.set_fontsize(8)
        for x in range(len(v.subset_labels)):
            if v.subset_labels[x] is not None:
                v.subset_labels[x].set_fontsize(8)
        plt.savefig(O2,dpi=300)
        plt.close(figV)

if __name__ == '__main__':
    main()
