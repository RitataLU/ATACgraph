#! /usr/bin/env python
#version 2 2019.9.6
import pysam
import sys
import subprocess
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy as sp
import scipy.fftpack
import numpy as np
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_bam')
    parser.add_argument('fragment_distribution_outname')
    parser.add_argument('fragment_fft_outname')
    args = parser.parse_args()
    input_bam = args.input_bam
    output1 = args.fragment_distribution_outname
    output2 = args.fragment_fft_outname
    fragment_distribution(input_bam,output1,output2)

def fragment_distribution(input_bam,output1,output2):
    #read bam file
    bam = pysam.AlignmentFile(input_bam)
    #get fragment lengths
    fragment_len = [ b.template_length for b in bam if b.template_length > 0 and
            b.template_length < 2000]

    #making histogram
    plt.style.use('classic')
    fig=plt.figure()
    ax1=fig.add_subplot(1,1,1)
    n, bins, patches = ax1.hist(fragment_len,bins=100,color='orange',range=(0,1000))
    ax1.xaxis.set_ticks_position('bottom')
    ax1.yaxis.set_ticks_position('left')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.xlabel('Fragment length', fontsize=20)
    plt.ylabel('Fragment count', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    ax1.set_title('Fragment distribution', fontsize=20)
    outfile1 = output1
    plt.savefig(outfile1, dpi=300,bbox_inches='tight')
    plt.close(fig)

    #calculating fft
    d2a = ffttable(n[15:96])
    d2b = ffttable(n[15:95])
    d2c = ffttable(n[15:94])
    d2d = ffttable(n[15:93])
    d2e = ffttable(n[15:92])
    d2f = ffttable(n[15:91])
    d2g = ffttable(n[15:90])
    d2h = ffttable(n[15:89])
    d2i = ffttable(n[15:88])
    d2j = ffttable(n[15:87])
    d2 = pd.concat([d2a,d2b,d2c,d2d,d2e,d2f,d2g,d2h,d2i,d2j],ignore_index=True)
    d2 = d2.sort_values(by=['freq'],ascending=False)
    max_y_pos = 1 / d2.loc[d2.value.idxmax(),'freq']
    plt.style.use('classic')
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(1/d2.freq,10 * np.log10(d2.value + 1))
    ax.axvline(max_y_pos,color = 'r',linestyle= '--')
    ax.set_xlim(0, 400)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    ax.set_xlabel('Period (bp)', fontsize=20)
    ax.set_ylabel('Power', fontsize=20)
    ax.set_title('Period of fragment distribution', fontsize=20)
    outfile2 = output2
    plt.savefig(outfile2,dpi=300,bbox_inches='tight')

def ffttable(selected):
    selected = np.log2(selected + 1)
    selected2 = [y-x for x, y in zip(selected, selected[1:])]
    fragment_fft = sp.fftpack.fft(selected2)
    fragment_psd = np.abs(fragment_fft) ** 2
    fftfreq = sp.fftpack.fftfreq(len(fragment_psd), 10.0)
    d = pd.DataFrame({'freq':fftfreq, 'value':fragment_psd})
    d2 = d[d.freq >0]
    return d2

if __name__ == '__main__':
    main()

