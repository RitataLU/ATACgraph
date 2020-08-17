#!/usr/bin/env python
#2019.12.9 added FC and pvalue 
#2929.03.09 added parameter for FC
import pybedtools
import argparse
import csv
import warnings
import pyBigWig
import numpy as np
from scipy.stats import ttest_ind

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_peakAs',help='Peaks files, seprate by comma')
    parser.add_argument('input_peakBs',help='Peaks files, seprate by comma')
    parser.add_argument('input_bwAs',help='BigWig files, seprate by comma')
    parser.add_argument('input_bwBs',help='BigWig files, seprate by comma')
    parser.add_argument('output_bedA')
    parser.add_argument('output_bedB')
    parser.add_argument('-c', '--fold_change', type=float, default=2.0,
        help='Fold change cutoff. Default:2')
    parser.add_argument('-p', '--p_value', type=float, default=0.05,
        help='P-value cutoff. Default:0.05')
    args = parser.parse_args()
#    inNameA = args.input_nameAs
#    inNameB = args.input_nameBs
    inPeakA = args.input_peakAs
    inPeakB = args.input_peakBs
    inBwA = args.input_bwAs
    inBwB = args.input_bwBs
    fold_change = args.fold_change
    p_value = args.p_value

#    inNameA = 'hESC_Ctrl_1_rmMul,hESC_Ctrl_2_rmMul'
#    inNameB = 'hESC_TKO_1_rmMul,hESC_TKO_2_rmMul'
#    groupNameA = inNameA.split(',')
#    groupNameB = inNameB.split(',')
    groupPeakA = inPeakA.split(',')
    groupPeakB = inPeakB.split(',')
    groupBwA = inBwA.split(',')
    groupBwB = inBwB.split(',')
#    groupPeakA = [name + '.peak_peaks.narrowPeak' for name in groupNameA]
#    groupPeakB = [name + '.peak_peaks.narrowPeak' for name in groupNameB]
#    groupBwA = [name + '_coverage.bw' for name in groupNameA]
#    groupBwB = [name + '_coverage.bw' for name in groupNameB]
    groupBw = groupBwA + groupBwB
    outBedA = args.output_bedA
    outBedB = args.output_bedB
    lenPeakA = len(groupPeakA)
    lenPeakB = len(groupPeakB)
    x = pybedtools.BedTool()
    result = x.multi_intersect(i=groupPeakA+groupPeakB,cluster='T')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore",category=FutureWarning)
        df = result.to_dataframe()

    df = df[df.name >= min(lenPeakA,lenPeakB)]
    for bwFile in groupBw:
        bw = pyBigWig.open(bwFile)
        df[bwFile] = [np.mean(bw.values(peak['chrom'],peak['start'] - 1,peak['end'])) for index,peak in df.iterrows()]

    df = df[df[groupBw].min(axis = 1) > 1]
    df.loc[:,'fc'] = np.log2(df[groupBwA].mean(axis=1)/df[groupBwB].mean(axis=1))
    if min(lenPeakA,lenPeakB) < 2:
        df.loc[:,'pvalue'] = 0
    else:
        df.loc[:,'pvalue'] = [ttest_ind(peak[groupBwA],peak[groupBwB],equal_var=False).pvalue for index,peak in df.iterrows()] 

    peaksA = df[(df.fc >= np.log2(fold_change)) & (df.pvalue < p_value)].copy()
    peaksA.loc[:,'id'] = ['PeakA_'+str(i) for i in range(1,len(peaksA)+1)]
    peaksA.loc[:,'fc2'] = ["%.2f" % var for var in peaksA.fc]
    peaksA = peaksA[['chrom','start','end','id','fc2']]
    peaksB = df[(df.fc <= -np.log2(fold_change)) & (df.pvalue < p_value)].copy()
    peaksB = df[(df.fc <= -1) & (df.pvalue < 0.05)].copy()
    peaksB.loc[:,'id'] = ['PeakB_'+str(i) for i in range(1,len(peaksB)+1)]
    peaksB.loc[:,'fc2'] = ["%.2f" % var for var in peaksB.fc]
    peaksB = peaksB[['chrom','start','end','id','fc2']]
    peaksA.to_csv(outBedA,sep="\t", index=False, header=None, quoting=csv.QUOTE_NONE)
    peaksB.to_csv(outBedB,sep="\t", index=False, header=None, quoting=csv.QUOTE_NONE)

if __name__ == '__main__':
    main()

