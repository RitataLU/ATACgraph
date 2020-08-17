### Integration multiomics data


import pandas as pd
import scipy.stats as ss

import subprocess, sys
import time
import argparse
import subprocess
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib_venn import venn2
from scipy.stats import hypergeom



def main():
        parser = argparse.ArgumentParser()
        parser.add_argument('atac_peak',help='ATAC-seq peak bed',type=str)
        parser.add_argument('other_peak',help='Other seq peak bed',type=str)
        parser.add_argument('ATAC_name',help='Name for ATAC peak',type=str)
        parser.add_argument('otherPeak_name',help='Name for other peak',type=str)
        parser.add_argument('overlap_name',help='Name for overlapping peak',type=str)
        parser.add_argument('Genome_size',help='Genome size(bp)',type=int)
        parser.add_argument('genes',help='Gene or promoter annotation file with bed6.bed format',type=str)
        

        args = parser.parse_args()

        A =args.atac_peak
        D =args.other_peak
        Gene = args.genes
        overPeak=args.overlap_name

        # A ='hESC_TKO_1_peaks.narrowPeak'
        # D = '/work1/home/yenmr/project/human_trim28/ChIP_public/KAP/ChIP_KAP_primed_peaks.narrowPeak'
        
        # Gene = 'hg19_gene_promoter_bed6.bed'
        #overPeak = 'TKO1_ChIP'
        # Genome_size = 3*10^9
        #pr53337opotion=args.propotion

        
        
        # integration peak finding (atac & chip )
        subprocess.call('''bedtools intersect -a %s -b %s -wo > %s'''%(A,D,overPeak+'.bed'), shell=True)

        #overlap regions
        subprocess.call('''bedtools intersect -a %s -b %s  > %s'''%(A,D,overPeak+'_region.bed'), shell=True)
        

        ## list of overlapping regions associate with genes 
        subprocess.call('''bedtools intersect -wo -a %s -b %s -wo | awk '{print$4}'> %s'''%(Gene,overPeak+'_region.bed', 'Integration_peaks_genes_list.txt'), shell=True)

        # bed name for integration peak
        o=str(overPeak)+'.bed'

        atac_peak=pd.read_csv(A,sep="\t",header=None)
        DNase_peak=pd.read_csv(D,sep="\t",header=None)
        oPeak=pd.read_csv(o,sep="\t",header=None)

        #calculate hypergeometric test
        ##mean: size mean
        Amean = (atac_peak[2]- atac_peak[1]).mean()
        Dmean = (DNase_peak[2]- DNase_peak[1]).mean()
        peaksize = (Amean+Dmean)/2
        
        #estimate peaks for whole genome
        g = int(args.Genome_size / peaksize)
        #g= int(Genome_size / peaksize) 
       
        # number of atac peaks
        a=len(atac_peak)
        # number of DNase peaks
        b=len(DNase_peak)
        # number of overlap peaks
        c=len(oPeak)
        
        hpd = ss.hypergeom(g, a, b)
        p = hpd.pmf(c)
#        print(format(p, ".3f"))


        otherPeak_name=args.otherPeak_name
        ATAC_name=args.ATAC_name

        out = venn2(subsets = (a, b, c), set_labels = (ATAC_name,otherPeak_name),set_colors=('purple', 'skyblue'), alpha = 0.7)
        
        if p < 1e-16:
                plt.text(0,0.6, 'p < 1e-16')
        else :
                plt.text(0,0.6, 'p =' )
                plt.text (0.15,0.6,p )

        
        for text in out.set_labels:
                text.set_fontsize(14)


        #save file
        mpl.rcParams['pdf.fonttype'] = 42
        mpl.rcParams['ps.fonttype'] = 42
        plt.savefig('Venndiagram_'+overPeak+'.pdf',dpi=300)

        subprocess.call('''rm *_region.bed*''', shell = True)

if __name__ == '__main__':
    main()
