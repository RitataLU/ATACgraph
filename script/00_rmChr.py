#!/usr/bin/env python
#from argparse import ArgumentParser

import pysam
import sys


def processInput (target,chrs):
    split = target.split(',') 
    a = []
    for eachTarget in split:
        if eachTarget in chrs:
            a = [x for x in chrs if x not in split]
        else:
            print "Chromosome {} is not in the list. Please check!".format(eachTarget)
            sys.exit()
    return (a,split)


def openBamFile(infile):
    try:
        samfile = pysam.AlignmentFile(infile,'rb')
    except IOError:
        sys.exit("The file {} does not exist".format(infile))
    except:
        sys.exit("The file {} does not have BAM format".format(infile))

    try:
#        samfile.has_index()
        samfile.check_index()
    except:
        sys.exit("{} does not indexed. You MUST index the file first!".format(infile))

    return samfile


def rmChr (infile,outfile,targetChr):
    samfile = openBamFile(infile)
    samout = pysam.AlignmentFile(outfile, "wb", template=samfile)
    chrs = samfile.references
    chrs_keep,chrs_rm = processInput (targetChr,chrs)
    keepread = 0
    for chrom in chrs_keep:
        fetchread = samfile.fetch(chrom)
        for read in fetchread:
            samout.write(read)
            keepread += 1
            if read.is_duplicate:
                dupread += 1

    totalread = sum([x.total for x in samfile.get_index_statistics()])
    removeread = totalread - keepread
    removeratio = float(removeread) / totalread

    for chrom in samfile.get_index_statistics():
        if chrom.contig in chrs_rm:
            rmcount = chrom.total
            print "Remove %s %d reads" % (chrom.contig,rmcount)

    print "Remove total %d out of %d (%.3f)" % (removeread,totalread,removeratio) 

    samout.close()
    samfile.close()


def main():
    if len(sys.argv) != 4:
        print "Usage: python 00_rmChr <input.bam> <output.bam> <chrM>"
        print ""
        print "If you need to remove multiple chromosomes, use comma"
        print "to separate. For example chrM,chrPt"
        sys.exit()

#    infile = 'Ctrl_1.bam'
#    outfile = "test.bam"
#    targetChr = 'chr1,chrX,chrY,chr20'
    infile = sys.argv[1]
    outfile = sys.argv[2]
    targetChr = sys.argv[3]
    rmChr(infile,outfile,targetChr)


if __name__ == '__main__':
    main()




