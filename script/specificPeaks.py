#!/usr/bin/env python
#2019.12.9 added FC and pvalue 
import pybedtools
import argparse
import csv
import warnings
import pyBigWig
import numpy as np
from scipy.stats import ttest_ind

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_nameAs',help='Sample names, seprate by comma')
    parser.add_argument('input_nameBs',help='Sample names, seprate by comma')
    parser.add_argument('output_bedA')
    parser.add_argument('output_bedB')
    args = parser.parse_args()
    inNameA = args.input_nameAs
    inNameB = args.input_nameBs
    groupNameA = inNameA.split(',')
    groupNameB = inNameB.split(',')
    groupPeakA = [name + '.peak_peaks.narrowPeak' for name in groupNameA]
    groupPeakB = [name + '.peak_peaks.narrowPeak' for name in groupNameB]
    groupBwA = [name + '_coverage.bw' for name in groupNameA]
    groupBwB = [name + '_coverage.bw' for name in groupNameB]
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

    if min(lenPeakA,lenPeakB) < 2:
        specA = ','.join(str(x) for x in range(1,lenPeakA +1))
        specB = ','.join(str(x) for x in range(lenPeakA +1,lenPeakA + lenPeakB +1))
        peaksA = df[df.score == specA].iloc[:,0:3]
        peaksB = df[df.score == specB].iloc[:,0:3]

    else:
        df = df[df.name > 1]
        for bwFile in groupBw:
            bw = pyBigWig.open(bwFile)
            df[bwFile] = [np.mean(bw.values(peak['chrom'],peak['start'] - 1,peak['end'])) for index,peak in df.iterrows()]
#        for bwFile in groupBwA:
#            bw = pyBigWig.open(bwFile)
#            df[bwFile] = [np.mean(bw.values(peak['chrom'],peak['start'] - 1,peak['end'])) for index,peak in df.iterrows()]

 #       for bwFile in groupBwB:
#            bw = pyBigWig.open(bwFile)
#            df[bwFile] = [np.mean(bw.values(peak['chrom'],peak['start'] - 1,peak['end'])) for index,peak in df.iterrows()]

        df = df[df[groupBw].min(axis = 1) > 1]
        df['fc'] = np.log2(df[groupBwA].mean(axis=1)/df[groupBwB].mean(axis=1))
        df['pvalue'] = [ttest_ind(peak[groupBwA],peak[groupBwB],equal_var=False).pvalue for index,peak in df.iterrows()] 
        peaksA = df[(df.fc >= 1) & (df.pvalue < 0.05)].iloc[:,0:3]
        peaksB = df[(df.fc <= -1) & (df.pvalue < 0.05)].iloc[:,0:3]
#        peaksA = df[(df.fc >= 1) & (df.pvalue < 0.05)]
#        peaksB = df[(df.fc <= -1) & (df.pvalue < 0.05)]

    peaksA.to_csv(outBedA,sep="\t", index=False, header=None, quoting=csv.QUOTE_NONE)
    peaksB.to_csv(outBedB,sep="\t", index=False, header=None, quoting=csv.QUOTE_NONE)

if __name__ == '__main__':
    main()

