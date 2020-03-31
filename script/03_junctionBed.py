#!/usr/bin/env python

import pysam
import pandas as pd
import numpy as np
import sys
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_bam')
    parser.add_argument('output_bed')
    parser.add_argument('-s','--separate',type = int,default=150, help="border length of long and short fragment(bp), default: 150 ")
    parser.add_argument('-b','--bin',type = int,default=10, help="bin size of a group on genome(bp), default: 10")
    parser.add_argument('-f','--filter',type = int,default=2, help="minimal  fragment junction tracks, default:2")
    args = parser.parse_args()
    input_bam = args.input_bam
    outbed = args.output_bed
    selectSize = args.separate
    filterFunc = args.filter
    bin_size = args.bin
    leftshift = 0
    rightshift = 0
    with open(outbed, 'w') as c:
        track_name = '''track name=junctions description="TopHat junctions"'''
        blank=""
        c.write("%s\n%s"%(track_name,blank))

    bam = pysam.AlignmentFile(input_bam)
    a = [[bam.references[b.rname],b.pos + leftshift, b.pos + b.tlen - rightshift] for b in bam if b.tlen > 0 and b.tlen < 2000]
    #chromosome name bam.references[b.rname]
    c = pd.DataFrame.from_records(a)
    c.columns = ['chr','str','end']
    c['str'] = c['str']//bin_size*bin_size
    c['end'] = c['end']//bin_size*bin_size
    c1_score = c.groupby(['chr','str','end']).size().reset_index()
    c1_score.columns = ['chr','str','end','count2']
    c1_score = c1_score[c1_score.count2 >= filterFunc]
    c1_score['location']=c1_score.end-c1_score.str-1
#    c1_score['size'] = c1_score['end'] - c1_score['str']
    c1_str = c1_score['str'].copy().apply(str)
    c1_end = c1_score['end'].copy().apply(str)
    c1_score['name'] = c1_score['chr'].str.cat(c1_str,sep ='_').str.cat(c1_end,sep = '_')
#    c1_score['dir']='+'
    c1_score['dir']=np.where(c1_score['location'] > selectSize,'+','-')
    c1_score['thickstart']=c1_score.str
    c1_score['thickend']=c1_score.end
    c1_score['block_count']=2
    c1_score['block_size']="1,1"
    c1_score['rgb']=np.where(c1_score['location'] > selectSize,'255,0,0','0,0,255')
    c1_score['zero']=0
    c1_score['block_location']=c1_score['zero'].map(str)+","+c1_score['location'].map(str)
    c1_junction_pd = c1_score.ix[:,('chr','str','end','name','count2','dir','thickstart','thickend','rgb','block_count','block_size','block_location')]
    c1_junction_pd.to_csv(outbed,mode='a', header=None, index=None, sep="\t")

if __name__ == '__main__':
    main()

