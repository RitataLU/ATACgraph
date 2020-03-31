#!/usr/bin/env python
#from argparse import ArgumentParser

import pysam
import sys
import argparse

def openBamFile(infile):
    try:
        infile = pysam.AlignmentFile(infile,'rb')
    except IOError:
        sys.exit("The file {} does not exist".format(infile))
    except:
        sys.exit("The file {} does not have BAM format".format(infile))

    try:
#        samfile.has_index()
        infile.check_index()
    except:
        sys.exit("{} does not indexed. You MUST index the file first!".format(infile))

    return infile


def selectRead (infile,outfile,selectMode,selectSize):
    samfile = openBamFile(infile)
#    bam = pysam.AlignmentFile(infile)
    samout = pysam.AlignmentFile(outfile, "wb", template=samfile)
    keepread = 0
    if selectMode == 1:
        for read in samfile:
            if abs(read.template_length) < selectSize:
                samout.write(read)
    else:
        for read in samfile:
            if abs(read.template_length) > selectSize:
                samout.write(read)

    samout.close()
    samfile.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_bam')
    parser.add_argument('output_bam')
    parser.add_argument('-f','--filter',type=int,default=150,help='fragment length selectione, dfault=150')
    parser.add_argument('-m','--mode',type=int,default=1,help='Select fragments smaller [1] or larger [2] than filter size. Default=1')
    args = parser.parse_args()
    infile = args.input_bam
    outfile = args.output_bam
    selectSize = args.filter
    selectMode = args.mode

    selectRead(infile,outfile,selectMode,selectSize)


if __name__ == '__main__':
    main()




