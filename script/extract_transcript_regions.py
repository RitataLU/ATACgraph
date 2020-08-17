#!/usr/bin/env python

# Stephen N. Floor
# 7 October 2014
# floor@berkeley.edu

import sys, os, argparse
from collections import defaultdict
import gzip
import re

class Transcript: 
    def __init__(self):
        #properties defined in UCSC knowngenes 
        self.name = '' 
        self.chrom = ''
        self.strand = '' 
        self.txStart = 0
        self.txEnd = 0 
        self.cdsStart = 0
        self.cdsEnd = 0
        self.exonCt = 0
        self.exonStarts = []
        self.exonEnds = []
        self.exonLengths = []
        
        #meta properties to be computed during construction.  these are lists of BED first four field tuples with the exception of Len terms which are the length of the total region for the gene 
        self.utr5 = []
        self.utr5Len = 0
        self.utr5start = []
        self.utr5startLen = 0
        self.cds = []
        self.cdsLen = 0 
        self.utr3 = []
        self.utr3Len = 0
        self.exons = []
        self.exonsLen = 0
        self.introns = []
        self.intronsLen = 0
    
        self.coding = False


    def __str__(self):  #currently roughly knownGenes format with a second line containing metadata 
        return "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n%s\t%d\t%s\t%d\t%s\t%d\t%s\t%d\t%s\t%d\t%s" % (self.name, self.chrom, self.strand, self.txStart, self.txEnd, self.cdsStart, self.cdsEnd, self.exonCt, self.exonStarts, self.exonEnds, self.utr5, self.utr5Len, self.cds, self.cdsLen, self.utr3, self.utr3Len, self.exons, self.exonsLen, self.introns, self.intronsLen, self.coding)
    
#BED format output is goal.  Fields are optional after featureEnd 
# chrom    featureStart   featureEnd   nameOfLine   score(0-1000)   strand   thickStart  thickEnd  itemRGBtuple  blockCount  blockSizes   blockStarts 

#this function returns a list of BED-formatted strings for the feature passed as region with multiple entries per region possible, one for each primitive (exon/intron) 
    def bedFormat(self, region="exons"):
        if (not self.coding and (region == "5utr" or region == "cds" or region == "3utr")):
            print ("Transcript.py bedFormat error: noncoding transcripts do not have 5utr/cds/3utr")
            return []

        returnVal = []

        if (region == "5utr"):
            for chunk in self.utr5:
                returnVal.append("%s\t%d\t%d\t%s\t0\t%c\t%d\t%d\t0" % (chunk[0], chunk[1], chunk[2], chunk[3]+"_5utr",self.strand, chunk[1],chunk[2]))

        elif (region == "5utr_start"):
            for chunk in self.utr5start:
                returnVal.append("%s\t%d\t%d\t%s\t0\t%c\t%d\t%d\t0" % (chunk[0], chunk[1], chunk[2], chunk[3]+"_5utr_start",self.strand, chunk[1],chunk[2]))

        elif (region == "cds"):
            for chunk in self.cds:
                returnVal.append("%s\t%d\t%d\t%s\t0\t%c\t%d\t%d\t0" % (chunk[0], chunk[1], chunk[2], chunk[3]+"_cds",self.strand, chunk[1],chunk[2]))

        elif (region == "3utr"):
            for chunk in self.utr3:
                returnVal.append("%s\t%d\t%d\t%s\t0\t%c\t%d\t%d\t0" % (chunk[0], chunk[1], chunk[2], chunk[3]+"_3utr",self.strand, chunk[1],chunk[2]))

        elif (region == "exons"):
            for chunk in self.exons:
                returnVal.append("%s\t%d\t%d\t%s\t0\t%c\t%d\t%d\t0" % (chunk[0], chunk[1], chunk[2], chunk[3]+"_exon",self.strand, chunk[1],chunk[2]))

        elif (region == "introns"):
            for chunk in self.introns:
                returnVal.append("%s\t%d\t%d\t%s\t0\t%c\t%d\t%d\t0" % (chunk[0], chunk[1], chunk[2], chunk[3]+"_intron",self.strand, chunk[1],chunk[2]))
        else:
            print ("Transcript.py bedFormat error: currently only regions 5utr/cds/3utr/exons/introns are supported")
            

        return returnVal

#BED format output is goal.  Fields are optional after featureEnd 
# chrom    featureStart   featureEnd   nameOfLine   score(0-1000)   strand   thickStart  thickEnd  itemRGBtuple  blockCount  blockSizes   blockStarts 

#blockCount - The number of blocks (exons) in the BED line.
#blockSizes - A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
#blockStarts - A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount. 

#this function returns a BED-formatted string for the feature passed as region with blocks defining the exons as per the BED file format 
    def blockBedFormat(self, region="exons"):
        if (not self.coding and (region == "5utr" or region == "cds" or region == "3utr")):
            print ("Transcript.py blockBedFormat error: noncoding transcripts do not have 5utr/cds/3utr")
            return ""

        returnVal = ""
        score = 0 
        rgb = 0

        if (region == "5utr"):
            chromStart = self.utr5[0][1] # start of feature is start of first block
            chromEnd = self.utr5[-1][2] # end of feature is end of last block 
            regionName = self.name + "_5utr"

            blockCount = len(self.utr5)
            blockSizes = ''.join(["%d," % (chunk[2]-chunk[1]) for chunk in self.utr5])
            blockStarts = ''.join(["%d," % (chunk[1]-chromStart) for chunk in self.utr5])

            #print "blockCount %d blockSizes %s blockStarts %s" % (blockCount, blockSizes,  blockStarts)

        elif (region == "5utr_start"):
            chromStart = self.utr5start[0][1] # start of feature is start of first block
            chromEnd = self.utr5start[-1][2] # end of feature is end of last block 
            regionName = self.name + "_5utr_start"

            blockCount = len(self.utr5start)
            blockSizes = ''.join(["%d," % (chunk[2]-chunk[1]) for chunk in self.utr5start])
            blockStarts = ''.join(["%d," % (chunk[1]-chromStart) for chunk in self.utr5start])

            #print "blockCount %d blockSizes %s blockStarts %s" % (blockCount, blockSizes,  blockStarts)
            
        elif (region == "cds"):
            chromStart = self.cds[0][1] # start of feature is start of first block
            chromEnd = self.cds[-1][2] # end of feature is end of last block 
            regionName = self.name + "_cds"
            blockCount = len(self.cds)
            blockSizes = ''.join(["%d," % (chunk[2]-chunk[1]) for chunk in self.cds])
            blockStarts = ''.join(["%d," % (chunk[1]-chromStart) for chunk in self.cds])

            #print "blockCount %d blockSizes %s blockStarts %s" % (blockCount, blockSizes,  blockStarts)

        elif (region == "3utr"):
            chromStart = self.utr3[0][1] # start of feature is start of first block
            chromEnd = self.utr3[-1][2] # end of feature is end of last block 
            regionName = self.name + "_3utr"
            blockCount = len(self.utr3)
            blockSizes = ''.join(["%d," % (chunk[2]-chunk[1]) for chunk in self.utr3])
            blockStarts = ''.join(["%d," % (chunk[1]-chromStart) for chunk in self.utr3])

            #print "blockCount %d blockSizes %s blockStarts %s" % (blockCount, blockSizes,  blockStarts)

        elif (region == "exons"):
            chromStart = self.exons[0][1] # start of feature is start of first block
            chromEnd = self.exons[-1][2] # end of feature is end of last block 
            regionName = self.name + "_exon"
            blockCount = len(self.exons)
            blockSizes = ''.join(["%d," % (chunk[2]-chunk[1]) for chunk in self.exons])
            blockStarts = ''.join(["%d," % (chunk[1]-chromStart) for chunk in self.exons])

            #print "blockCount %d blockSizes %s blockStarts %s" % (blockCount, blockSizes,  blockStarts)

        elif (region == "introns"):
            chromStart = self.introns[0][1] # start of feature is start of first block
            chromEnd = self.introns[-1][2] # end of feature is end of last block 
            regionName = self.name + "_intron"
            blockCount = len(self.introns)
            blockSizes = ''.join(["%d," % (chunk[2]-chunk[1]) for chunk in self.introns])
            blockStarts = ''.join(["%d," % (chunk[1]-chromStart) for chunk in self.introns])

            #print "blockCount %d blockSizes %s blockStarts %s" % (blockCount, blockSizes,  blockStarts)
        else:
            print ("UCSCKnownGene blockBedFormat error: currently only regions 5utr/cds/3utr/exons/introns are supported")
            
        returnVal = "%s\t%d\t%d\t%s\t%d\t%c\t%d\t%d\t%s\t%d\t%s\t%s" % (self.chrom, chromStart, chromEnd, regionName, score, self.strand, chromStart, chromEnd, rgb, blockCount, blockSizes, blockStarts)

        return returnVal
        
    def computeMetadata(self): 
    # -- begin computing metadata -- 

    # -- note: chose clarity of code and conditionals here over most efficient computation (i.e. some clauses may be redundant)

        if (self.strand == "+"): 
        #print ("DBUG - exonCt %d i %d exonEnds[i] %d cdsStart %d exonStarts[i] %d cdsEnd %d") % \
            #    (self.exonCt, i, self.exonEnds[i], self.cdsStart, self.exonStarts[i], self.cdsEnd)
            for i in range (self.exonCt): 
                if (self.cdsStart != self.cdsEnd): # if this is a coding transcript
                    self.coding = True
                # -- first compute 5'utr, CDS, 3'utr regions --
                #case 1 - exon spans 5' UTR/CDS/3' UTR
                    if (self.exonStarts[i] < self.cdsStart and self.exonEnds[i] > self.cdsEnd):
                        self.utr5.append((self.chrom, self.exonStarts[i], self.cdsStart, self.name))
                        self.utr5Len += self.cdsStart - self.exonStarts[i]
                        self.utr5start.append((self.chrom, self.exonStarts[i], self.cdsStart, self.name)) # for now just append the 5' utr exons to the utr5start 
                        self.utr5startLen += self.cdsStart - self.exonStarts[i]
                        self.cds.append((self.chrom, self.cdsStart, self.cdsEnd, self.name))
                        self.cdsLen += self.cdsEnd - self.cdsStart
                        self.utr3.append((self.chrom, self.cdsEnd, self.exonEnds[i], self.name))
                        self.utr3Len += self.exonEnds[i] - self.cdsEnd
                #case 2 - exon spans 5' UTR/CDS junction
                    elif (self.exonStarts[i] < self.cdsStart and self.exonEnds[i] >= self.cdsStart):
                        self.utr5.append((self.chrom, self.exonStarts[i], self.cdsStart, self.name))
                        self.utr5Len += self.cdsStart - self.exonStarts[i]
                        self.utr5start.append((self.chrom, self.exonStarts[i], self.cdsStart, self.name)) 
                        self.utr5startLen += self.cdsStart  - self.exonStarts[i]
                        self.cds.append((self.chrom, self.cdsStart, self.exonEnds[i], self.name))
                        self.cdsLen += self.exonEnds[i]- self.cdsStart
                #case 3 - exon spans CDS/3'UTR junction 
                    elif (self.exonStarts[i] >= self.cdsStart and self.exonStarts[i] <= self.cdsEnd and self.exonEnds[i] > self.cdsEnd):
                        self.cds.append((self.chrom, self.exonStarts[i], self.cdsEnd, self.name))
                        self.cdsLen += self.cdsEnd - self.exonStarts[i]
                        self.utr3.append((self.chrom, self.cdsEnd, self.exonEnds[i], self.name))
                        self.utr3Len += self.exonEnds[i] - self.cdsEnd
                #case 4 - exon is 5' UTR only 
                    elif (self.exonStarts[i] < self.cdsStart and self.exonEnds[i] < self.cdsStart): 
                        self.utr5.append((self.chrom, self.exonStarts[i], self.exonEnds[i], self.name))
                        self.utr5Len += self.exonEnds[i] - self.exonStarts[i]
                        self.utr5start.append((self.chrom, self.exonStarts[i], self.exonEnds[i], self.name)) 
                        self.utr5startLen += self.exonEnds[i] - self.exonStarts[i]
                #case 5 - exon is CDS only
                    elif (self.exonStarts[i] >= self.cdsStart and self.exonEnds[i] <= self.cdsEnd):
                        self.cds.append((self.chrom, self.exonStarts[i], self.exonEnds[i], self.name))
                        self.cdsLen += self.exonEnds[i] - self.exonStarts[i]
                #case 6 - exon is 3' UTR only 
                    elif (self.exonStarts[i] > self.cdsEnd and self.exonEnds[i] > self.cdsEnd):
                        self.utr3.append((self.chrom, self.exonStarts[i], self.exonEnds[i], self.name))
                        self.utr3Len += self.exonEnds[i] - self.exonStarts[i]
                    else: 
                        print ("Thar be dragons - Transcript computeMetadata + stranded gene region parsing")


            # -- generate combined exonic and intronic regions -- 
            #exons are easy 
                self.exons.append((self.chrom, self.exonStarts[i], self.exonEnds[i], self.name))
                self.exonsLen += self.exonEnds[i] - self.exonStarts[i]
            
            #print "DBUG2: i %d self.exonCt-1 %d self.exonEnds %s self.exonStarts %s" % (i, self.exonCt-1, self.exonEnds, self.exonStarts)
        
                if (i < self.exonCt - 1): # only compute introns for nonterminal exons
                # an intron is the region between the end of the current exon and start of the next 
                    self.introns.append((self.chrom, self.exonEnds[i], self.exonStarts[i+1], self.name))
                    self.intronsLen += self.exonStarts[i+1] - self.exonEnds[i] 

        # append 27 nt of the coding sequence onto the utr5start regions by creating another block, only if there is a 5' utr already 
            if (self.coding and self.utr5Len > 0): 
                added = 0 
                index = 0
            
                while(added < 27 and index < len(self.cds)):
                    if (self.cds[index][2] - self.cds[index][1]) > 27 - added: 
                        self.utr5start.append((self.chrom, self.cds[index][1], self.cds[index][1] + 27 - added, self.name))
                        self.utr5startLen += 27 - added
                        added += 27 - added
                    else: 
                        self.utr5start.append((self.chrom, self.cds[index][1], self.cds[index][2], self.name))
                        self.utr5startLen += self.cds[index][2] - self.cds[index][1]
                        added += self.cds[index][2] - self.cds[index][1]
                        index += 1

                if (added < 27):
                    #print "Transcript.py: aborting 5' UTR start for transcript %s because %d nts added (CDS length %d)" % (self.name, added, self.cdsLen) 
                    self.utr5startLen = 0
                    self.utr5start = []

        elif (self.strand == "-"):
     #uc001ach.2	    chr1    -	    910578  917473  911551  916546  5	    910578,911878,914260,916516,917444,	    911649,912004,916037,916553,917473,	    Q5SV97  uc001ach.2
            #	name		chrom	strand	txStart txEnd	cdsStart self.cdsEnd exonCt	exonStarts		exonEnds		proteinID  alignID 
            # for the minus strand everything is the same except the order of encountering regions is reversed
            # i.e. 3' UTR -> CDS -> 5' UTR 
            
            for i in range (self.exonCt): 
            #print ("DBUG - exonCt %d i %d self.exonEnds[i] %d self.cdsStart %d exonStarts[i] %d self.cdsEnd %d") % \
                #    (self.exonCt, i, self.exonEnds[i], self.cdsStart, self.exonStarts[i], self.cdsEnd)
                
                if (self.cdsStart != self.cdsEnd):
                    self.coding = True 
                # -- first compute 5'utr, CDS, 3'utr regions --
                # -- this is the same as for + sense except 5' UTR and 3' UTR are swapped throughout
                #case 1 - exon spans 3' UTR/CDS/5' UTR
                    if (self.exonStarts[i] < self.cdsStart and self.exonEnds[i] > self.cdsEnd):
                        self.utr3.append((self.chrom, self.exonStarts[i], self.cdsStart, self.name))
                        self.utr3Len += self.cdsStart - self.exonStarts[i]
                        self.cds.append((self.chrom, self.cdsStart, self.cdsEnd, self.name))
                        self.cdsLen += self.cdsEnd - self.cdsStart
                        self.utr5.append((self.chrom, self.cdsEnd, self.exonEnds[i], self.name))
                        self.utr5Len += self.exonEnds[i] - self.cdsEnd
                        self.utr5start.append((self.chrom, self.cdsEnd, self.exonEnds[i], self.name))
                        self.utr5startLen += self.exonEnds[i] - (self.cdsEnd)
                #case 2 - exon spans 3' UTR/CDS junction
                    elif (self.exonStarts[i] < self.cdsStart and self.exonEnds[i] >= self.cdsStart):
                        self.utr3.append((self.chrom, self.exonStarts[i], self.cdsStart, self.name))
                        self.utr3Len += self.cdsStart - self.exonStarts[i]
                        self.cds.append((self.chrom, self.cdsStart, self.exonEnds[i], self.name))
                        self.cdsLen += self.exonEnds[i]- self.cdsStart
                #case 3 - exon spans CDS/5'UTR junction 
                    elif (self.exonStarts[i] >= self.cdsStart and self.exonStarts[i] <= self.cdsEnd and self.exonEnds[i] > self.cdsEnd):
                        self.cds.append((self.chrom, self.exonStarts[i], self.cdsEnd, self.name))
                        self.cdsLen += self.cdsEnd - self.exonStarts[i]
                        self.utr5.append((self.chrom, self.cdsEnd, self.exonEnds[i], self.name))
                        self.utr5Len += self.exonEnds[i] - self.cdsEnd
                        self.utr5start.append((self.chrom, self.cdsEnd, self.exonEnds[i], self.name))
                        self.utr5startLen += self.exonEnds[i] - (self.cdsEnd)
                #case 4 - exon is 3' UTR only 
                    elif (self.exonStarts[i] < self.cdsStart and self.exonEnds[i] < self.cdsStart): 
                        self.utr3.append((self.chrom, self.exonStarts[i], self.exonEnds[i], self.name))
                        self.utr3Len += self.exonEnds[i] - self.exonStarts[i]
                #case 5 - exon is CDS only
                    elif (self.exonStarts[i] >= self.cdsStart and self.exonEnds[i] <= self.cdsEnd):
                        self.cds.append((self.chrom, self.exonStarts[i], self.exonEnds[i], self.name))
                        self.cdsLen += self.exonEnds[i] - self.exonStarts[i]
                #case 6 - exon is 5' UTR only 
                    elif (self.exonStarts[i] > self.cdsEnd and self.exonEnds[i] > self.cdsEnd):
                        self.utr5.append((self.chrom, self.exonStarts[i], self.exonEnds[i], self.name))
                        self.utr5Len += self.exonEnds[i] - self.exonStarts[i]
                        self.utr5start.append((self.chrom, self.exonStarts[i] , self.exonEnds[i], self.name))
                        self.utr5startLen += self.exonEnds[i] - self.exonStarts[i]
                    else: 
                        print ("Thar be dragons - Transcript computeMetadata - stranded gene region parsing")
                    
            #else: 
            #    print "- strand noncoding transcript"
                

            # -- generate combined exonic and intronic regions -- 
            #exons are easy 
                self.exons.append((self.chrom, self.exonStarts[i], self.exonEnds[i], self.name))
                self.exonsLen += self.exonEnds[i] - self.exonStarts[i]
            
                if (i < self.exonCt - 1): # only compute introns for nonterminal exons
                # an intron is the region between the end of the current exon and start of the next 
                    self.introns.append((self.chrom, self.exonEnds[i], self.exonStarts[i+1], self.name))
                    self.intronsLen += self.exonStarts[i+1] - self.exonEnds[i] 
                
        # append 27 nt of the coding sequence onto the utr5start regions by creating another block, only if there is a 5' utr already 
            if (self.coding and self.utr5Len > 0): 
                added = 0 
                index = -1
                
                while(added < 27 and index >= len(self.cds)*-1):  #cdsEnd is the start, and the last exon is the first exon of the cds 
                #print self.cds
                #print self.cds[index]
                
                # need to insert at the beginning here and not append
                    
                    if (self.cds[index][2] - self.cds[index][1]) > 27 - added: 
                        self.utr5start.insert(0,(self.chrom, self.cds[index][2] - (27 - added), self.cds[index][2], self.name))
                        self.utr5startLen += 27 - added 
                        added += 27 - added
                    else: 
                        self.utr5start.insert(0,(self.chrom, self.cds[index][1], self.cds[index][2], self.name))
                        self.utr5startLen += self.cds[index][2] - self.cds[index][1]
                        added += self.cds[index][2] - self.cds[index][1]
                        index -= 1

                if (added < 27):
                    #print "Transcript.py: aborting 5' UTR start for transcript %s because %d nts added (CDS length %d)" % (self.name, added, self.cdsLen) 
                    self.utr5startLen = 0
                    self.utr5start = []
                        
        else:
            print ("Thar be dragons - Transcript computeMetadata strand does not match + or -")
        

#example line format from knownGenes file (from UCSC) 
    #    uc010nxq.1      chr1    +       11873   14409   12189   13639   3       11873,12594,13402,      12227,12721,14409,      B7ZGX9  uc010nxq.1
    # line format
    #    name            chrom   strand  txStart txEnd   cdsStart cdsEnd exonCt  exonStarts              exonEnds                proteinID  alignID 

def createUCSCTranscript(knownGeneLine):
    foo = Transcript()
    line = knownGeneLine.split()

    # -- read in knownGene fields -- 

    foo.name = line[0]
    foo.chrom = line[1]
    foo.strand = line[2]
    foo.txStart = int(line[3])
    foo.txEnd = int(line[4])
    foo.cdsStart = int(line[5])
    foo.cdsEnd = int(line[6])
    foo.exonCt = int(line[7])
    
    starts = line[8].split(",")
    ends = line[9].split(",")

    for i in range(foo.exonCt):
        foo.exonStarts.append(int(starts[i]))
        foo.exonEnds.append(int(ends[i]))

    foo.computeMetadata() 
        
    return foo 

# input to createGTFTranscript below must be a list of dictionaries for each line of the input GTF file 
# these are created inside knowngenes_to_transcript_regions.py 

# example input: 

#[{'gene_name': 'DDX11L1', 'seqname': '1', 'end': '12227', 'start': '11869', 'frame': None, 'transcript_source': 'havana', 'feature': 'exon', 'exon_number': '1', 'exon_id': 'ENSE00002234944', 'tss_id': 'TSS15145', 'source': 'processed_transcript', 'gene_source': 'ensembl_havana', 'score': None, 'gene_biotype': 'pseudogene', 'gene_id': 'ENSG00000223972', 'transcript_id': 'ENST00000456328', 'transcript_name': 'DDX11L1-002', 'strand': '+'}, {'seqname': '1', 'end': '14409', 'start': '11869', 'frame': None, 'transcript_source': 'havana', 'feature': 'transcript', 'gene_id': 'ENSG00000223972', 'tss_id': 'TSS15145', 'source': 'processed_transcript', 'gene_source': 'ensembl_havana', 'score': None, 'gene_biotype': 'pseudogene', 'gene_name': 'DDX11L1', 'transcript_id': 'ENST00000456328', 'transcript_name': 'DDX11L1-002', 'strand': '+'}]

# keys for each dict:
#  gene_name
#  seqname
#  start
#  end
#  frame
#  transcript_source
#  feature
#  exon_number
#  exon_id
#  tss_id
#  source
#  gene_source
#  score
#  gene_biotype
#  gene_id
#  transcript_id
#  transcript_name
#  strand

def createGTFTranscript(gtfLines):
    foo = Transcript()


    # these properties (better be) all identical for each entry in the list of dicts 

    first = gtfLines[0] 

    foo.name = first["transcript_id"]
    foo.chrom = first["seqname"]
    foo.strand = first["strand"]

    # now process all lines for this transcript ID 

    for dict in gtfLines: 

        # ensembl GTFs have special lines where feature = "transcript" and feature = "CDS" that define the transcript and CDS start/ends, respectively 

        # GTF files are closed intervals while BED are right-open-left-closed, so --- 
        #   need to subtract one from all start coordinates? seems counterintuitive maybe the input genome.fa is zero based? 

        if (dict["feature"] == "transcript"):
            if (foo.txStart == 0 or int(dict["start"]) < foo.txStart):
                foo.txStart = int(dict["start"]) - 1
            if (foo.txEnd == 0 or int(dict["end"]) > foo.txEnd):
                foo.txEnd = int(dict["end"])


        if (dict["feature"] == "CDS"):
            if (foo.cdsStart == 0 or int(dict["start"]) < foo.cdsStart):
                foo.cdsStart = int(dict["start"]) - 1
            if (foo.cdsEnd== 0 or int(dict["end"]) > foo.cdsEnd):
                foo.cdsEnd = int(dict["end"]) 

        if (dict["feature"] == "exon"):
            foo.exonCt += 1

            foo.exonStarts.append(int(dict["start"]) - 1)
            foo.exonEnds.append(int(dict["end"]))

    foo.exonStarts = sorted(foo.exonStarts)
    foo.exonEnds = sorted(foo.exonEnds) 

    foo.computeMetadata() 



    return foo


#GTF
GTF_HEADER  = ['seqname', 'source', 'feature', 'start', 'end', 'score',
           'strand', 'frame']
R_SEMICOLON = re.compile(r'\s*;\s*')
R_COMMA     = re.compile(r'\s*,\s*')
R_KEYVALUE  = re.compile(r'(\s+|\s*=\s*)')

def lines(filename):
    """Open an optionally gzipped GTF file and generate a dict for each line.
    """
    fn_open = gzip.open if filename.endswith('.gz') else open

    with fn_open(filename) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            else:
                yield parse(line)


def parse(line):
    """Parse a single GTF line and return a dict.
    """
    result = {}

    fields = line.rstrip().split('\t')

    for i, col in enumerate(GTF_HEADER):
        result[col] = _get_value(fields[i])

    # INFO field consists of "key1=value;key2=value;...".
    infos = re.split(R_SEMICOLON, fields[8])

    for i, info in enumerate(infos, 1):
        # It should be key="value".
        try:
            key, _, value = re.split(R_KEYVALUE, info)
        # But sometimes it is just "value".
        except ValueError:
            key = 'INFO{}'.format(i)
            value = info
        # Ignore the field if there is no value.
        if value:
            result[key] = _get_value(value)

    return result


def _get_value(value):
    if not value:
        return None

    # Strip double and single quotes.
    value = value.strip('"\'')

    # Return a list if the value has a comma.
    if ',' in value:
        value = re.split(R_COMMA, value)
    # These values are equivalent to None.
    elif value in ['', '.', 'NA']:
        return None

    return value


#from Transcript import *

print (" ----------------------------------")
print ("| Extract Regions from annotations |")
print ("|  snf        Fall 2014            |")
print (" ----------------------------------\n\n")


# ------ ARGUMENT PARSING ----------

parser = argparse.ArgumentParser(description="Create transcript regions (5' UTR/CDS/3'UTR etc) from knownGenes or a GTF") 

parser.add_argument("-i", "--input", help="input filename", required=True)
parser.add_argument("-o", "--output", help="output basename", required=True) 
parser.add_argument("--ucsc", help="Read from a UCSC knownGenes formatted file (BED)", action="store_true")
parser.add_argument("--gtf", help="Read from a GTF (only tested with Ensembl GTFs)", action="store_true")

args = parser.parse_args() 

if ( (not (args.ucsc or args.gtf)) or (args.ucsc and args.gtf)):
    sys.exit("FATAL: must set one but not both of --ucsc and --gtf") 

# output filenames:
utr5FName = args.output + "_5utr.bed"
utr5StartFName = args.output + "_5utr_start.bed"
cdsFName = args.output + "_cds.bed"
utr3FName = args.output + "_3utr.bed"
exonFName = args.output + "_exons.bed"
intronFName = args.output + "_introns.bed"
codingExonFName = args.output + "_codingexons.bed"
codingIntronFName = args.output + "_codingintrons.bed" # note that these are introns from coding genes, not necessarily introns that make it to mRNA
noncodingExonFName = args.output + "_noncodingexons.bed"
noncodingIntronFName = args.output + "_noncodingintrons.bed"

#keep track of where we are
genesRead = 0

# parameters that should be passed via the cmd line
useBlocks = True


#terminate if output files exist

#if os.path.exists(utr5FName) or os.path.exists(utr5StartFName) or os.path.exists(cdsFName) or os.path.exists(utr3FName) or os.path.exists(exonFName) or os.path.exists(intronFName) \
#        or os.path.exists(codingExonFName) or os.path.exists(codingIntronFName) or os.path.exists(noncodingExonFName) or os.path.exists(noncodingIntronFName):
#    sys.exit("ERROR: output basename %s files already exist" % args.output)

#process the file

with open(utr5FName, "w") as utr5File, open(utr5StartFName, "w") as utr5StartFile, open(cdsFName, "w") as cdsFile, \
        open(utr3FName, "w") as utr3File, open(exonFName, "w") as exonFile, open (intronFName, "w") as intronFile, \
        open(codingExonFName, "w") as codingExonFile, open(codingIntronFName, "w") as codingIntronFile, \
        open(noncodingExonFName, "w") as noncodingExonFile, open(noncodingIntronFName, "w") as noncodingIntronFile:

    def writeOutput(gene):
        if (useBlocks): # output all region primitives on the same line by specifying nBlocks and lists inside the BED output
            if(gene.coding):
                #blockBedFormat is one line by definition
                if (gene.utr5Len > 0): utr5File.write(gene.blockBedFormat(region="5utr") + "\n")
                if (gene.utr5startLen > 0): utr5StartFile.write(gene.blockBedFormat(region="5utr_start") + "\n")
                if (gene.cdsLen > 0): cdsFile.write(gene.blockBedFormat(region="cds") + "\n")
                if (gene.utr3Len > 0): utr3File.write(gene.blockBedFormat(region="3utr") + "\n")
            
                if (gene.exonsLen > 0):
                    exonFile.write(gene.blockBedFormat(region="exons") + "\n")
                    codingExonFile.write(gene.blockBedFormat(region="exons") + "\n")
                
                if (gene.intronsLen > 0):
                    intronFile.write(gene.blockBedFormat(region="introns") + "\n")
                    codingIntronFile.write(gene.blockBedFormat(region="introns") + "\n")
                    
            else: # noncoding transcripts just have exons and introns
                if (gene.exonsLen > 0):
                    exonFile.write(gene.blockBedFormat(region="exons") + "\n")
                    noncodingExonFile.write(gene.blockBedFormat(region="exons") + "\n")

                if (gene.intronsLen > 0):
                    intronFile.write(gene.blockBedFormat(region="introns") + "\n")
                    noncodingIntronFile.write(gene.blockBedFormat(region="introns") + "\n")

        else: # output one line per region primitive instead of combining regions via blocks
            if(gene.coding):
                for entry in gene.bedFormat(region="5utr"):
                    utr5File.write(entry + "\n")
                for entry in gene.bedFormat(region="5utr_start"):
                    utr5StartFile.write(entry + "\n")
                for entry in gene.bedFormat(region="cds"):
                    cdsFile.write(entry + "\n")
                for entry in gene.bedFormat(region="3utr"):
                    utr3File.write(entry + "\n")

                for entry in gene.bedFormat(region="exons"):
                    exonFile.write(entry + "\n")
                    codingExonFile.write(entry + "\n")

                for entry in gene.bedFormat(region="introns"):
                    intronFile.write(entry + "\n")
                    codingIntronFile.write(entry + "\n")

            else: # noncoding transcripts just have exons and introns
                for entry in gene.bedFormat(region="exons"):
                    exonFile.write(entry + "\n")
                    noncodingExonFile.write(entry + "\n")

                for entry in gene.bedFormat(region="introns"):
                    intronFile.write(entry + "\n")
                    noncodingIntronFile.write(entry + "\n")


    if (args.ucsc): 
        with open(args.input, "r") as genesFile: 
        
            for line in genesFile:
                # all of the knowngenes parsing and metadata construction is done inside Transcript.py, especially the createGene method

                gene = createUCSCTranscript(line) 
                genesRead += 1

                writeOutput(gene)

                if (not genesRead % 2500):
                    print ("Processed %d entries..." %  genesRead)

                
    elif (args.gtf): 
            
            # first parse the entire file into a dictionary of lists

        txDict = defaultdict(list) 

        print ("Building GTF dictionary...")

        # the issue here is that lines for various transcripts may be interleaved, so can either create lots of objects, or a giant dict. opted for giant dict. 
        for line in lines(args.input): 
            # only want to read in lines corresponding to these features
            if line["feature"] in ["exon", "CDS", "start_codon", "stop_codon"]:
                txDict[line["transcript_id"]].append(line)
                genesRead += 1

                if (not genesRead % 25000):
                    print ("\tProcessed %d lines..." %  genesRead)

        print ("Dictionary built.")

        print ("Writing transcript properties.")
        genesRead = 0
        
        # now create a Transcript object for each transcript and output it 

        for key in txDict: 

            #print key

            tx = createGTFTranscript(txDict[key])

            #print tx 
            writeOutput(tx)
            genesRead += 1
            
            if (not genesRead % 2500):
                print ("\tProcessed %d entries..." %  genesRead)


print ("Processed %d entries." %  genesRead)




