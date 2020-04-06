#!/usr/bin/env python
#2920.03.16 
import pybedtools
import argparse
import os
import sys, subprocess
import csv
import warnings
import pyBigWig
import numpy as np
from scipy.stats import ttest_ind

ATACgraph_dir = os.path.dirname(os.path.abspath(__file__))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_peakAs',help='Peaks files, seprate by comma')
    parser.add_argument('input_peakBs',help='Peaks files, seprate by comma')
    parser.add_argument('output_bedA')
    parser.add_argument('output_bedB')
#    parser.add_argument('-p', '--p_value', type=float, default=0.05,
#        help='P-value cutoff. Default:0.05')
    args = parser.parse_args()
    inPeakA = args.input_peakAs
    inPeakB = args.input_peakBs
    inPeakA = 'hESC_Ctrl_1_rmMul_whole.peak_peaks.narrowPeak,hESC_Ctrl_2_rmMul_whole.peak_peaks.narrowPeak'
    inPeakB = 'hESC_TKO_1_rmMul_whole.peak_peaks.narrowPeak,hESC_TKO_2_rmMul_whole.peak_peaks.narrowPeak'
    p_value = args.p_value
    groupPeakA = inPeakA.split(',')
    groupPeakB = inPeakB.split(',')
    outBedA = args.output_bedA
    outBedB = args.output_bedB
    lenPeakA = len(groupPeakA)
    lenPeakB = len(groupPeakB)
    idrAFile = "IDR_A.txt"
    idrBFile = "IDR_B.txt"
    idrACmd = '''idr --sample %s --output-file %s'''%(" ".join(groupPeakA),idrAFile)
    idrBCmd = '''idr --sample %s --output-file %s'''%(" ".join(groupPeakB),idrBFile)
    specACmd = '''bedtools subtract -a %s -b %s -A > %s'''%(idrAFile,idrBFile,outBedA)
    specBCmd = '''bedtools subtract -a %s -b %s -A > %s'''%(idrBFile,idrAFile,outBedB)
    print(idrACmd)
    subprocess.call(idrACmd, shell=True)
    print(idrBCmd)
    subprocess.call(idrBCmd, shell=True)
    print(specACmd)
    subprocess.call(specACmd, shell=True)
    print(specBCmd)
    subprocess.call(specBCmd, shell=True)


if __name__ == '__main__':
    main()

