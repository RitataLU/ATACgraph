#!/usr/bin/env python
# coding: utf-8

import argparse
import warnings
import csv
import pyBigWig
import numpy as np
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from scipy.stats import norm, multivariate_normal
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_bed',help='Promoter bed file')
    parser.add_argument('input_bwAs',help='BigWig files, seprate by comma')
    parser.add_argument('input_bwBs',help='BigWig files, seprate by comma')
    parser.add_argument('output_bedA')
    parser.add_argument('output_bedB')
    parser.add_argument('-m', '--minimum_abundance', type=float, default=10.0,
        help='Minimum abundance of opened promoter. Default:10.0')
    parser.add_argument('-p', '--p_value', type=float, default=0.05,
        help='P-value cutoff. Default:0.05')
    args = parser.parse_args()
    minOpenAbundance = args.minimum_abundance
    prob_cutoff = args.p_value
    #inBwA = 'hESC_Ctrl_1_rmMul_whole_coverage.bw,hESC_Ctrl_2_rmMul_whole_coverage.bw'
    #inBwB = 'hESC_TKO_1_rmMul_whole_coverage.bw,hESC_TKO_2_rmMul_whole_coverage.bw'
    inBwA = args.input_bwAs
    inBwB = args.input_bwBs
    bed = args.input_bed
    groupBwA = inBwA.split(',')
    groupBwB = inBwB.split(',')
    groupBw = groupBwA + groupBwB
    outBedA = args.output_bedA
    outBedB = args.output_bedB
    df = pd.read_csv(bed,sep="\t",header=None,names=['chrom','start','end','id'],usecols=[0,1,2,3])
    df = df.dropna()
    #remove chromosome that are not in bw
    bw0 = pyBigWig.open(groupBw[0])
    chrs = bw0.chroms().keys()
    df = df[df.chrom.isin(chrs)]
    for bwFile in groupBw:
        bw = pyBigWig.open(bwFile)
        df[bwFile] = [np.mean(bw.values(gene['chrom'],gene['start'] - 1,gene['end'])) for index,gene in df.iterrows()]  

    log_prob_cutoff = np.log(prob_cutoff)
    df = df[(df[groupBw].min(axis = 1) > 0) & (df[groupBw].max(axis = 1) > minOpenAbundance)]
    df.loc[:,'fc'] = np.log2(df[groupBwA].mean(axis=1)/df[groupBwB].mean(axis=1))    
    X = df.loc[:,'fc'].values.reshape(-1, 1)
    componemts = 1
    gmm = GaussianMixture(n_components=componemts,covariance_type='full').fit(X)
    b=gmm.score_samples(X)
    geneA = df[(b < log_prob_cutoff) & (df.fc > 0)].copy()
    geneB = df[(b < log_prob_cutoff) & (df.fc < 0)].copy()
    geneA.loc[:,'fc2'] = ["%.2f" % var for var in geneA.fc]
    geneB.loc[:,'fc2'] = ["%.2f" % var for var in geneB.fc]
    geneA = geneA[['chrom','start','end','id','fc2']]
    geneB = geneB[['chrom','start','end','id','fc2']]
    geneA.to_csv(outBedA,sep="\t", index=False, header=None, quoting=csv.QUOTE_NONE)
    geneB.to_csv(outBedB,sep="\t", index=False, header=None, quoting=csv.QUOTE_NONE)



if __name__ == '__main__':
    main()

