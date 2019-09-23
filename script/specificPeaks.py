#!/usr/bin/env python

import pybedtools
import argparse
import csv
import warnings

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_peakAs',help='Peak beds, seprate by comma')
    parser.add_argument('input_peakBs',help='Peak beds, seprate by comma')
    parser.add_argument('output_bedA')
    parser.add_argument('output_bedB')
    args = parser.parse_args()
    inA = args.input_peakAs
    inB = args.input_peakBs
    outBedA = args.output_bedA
    outBedB = args.output_bedB
    groupA = inA.split(',')
    groupB = inB.split(',')
    lenA = len(groupA)
    lenB = len(groupB)
    x = pybedtools.BedTool()
    result = x.multi_intersect(i=groupA+groupB,cluster='T')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore",category=FutureWarning)
        df = result.to_dataframe()

#    df = result.to_dataframe()
    specA = ','.join(str(x) for x in range(1,lenA +1))
    specB = ','.join(str(x) for x in range(lenA +1,lenA + lenB +1))
    peaksA = df[df.score == specA].iloc[:,0:3]
    peaksB = df[df.score == specB].iloc[:,0:3]

    peaksA.to_csv(outBedA,sep="\t", index=False, header=None, quoting=csv.QUOTE_NONE)
    peaksB.to_csv(outBedB,sep="\t", index=False, header=None, quoting=csv.QUOTE_NONE)

if __name__ == '__main__':
    main()

