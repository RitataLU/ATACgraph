
# callpeak & make bw file seperate only integration, no integration
import sys
import os
import subprocess, sys
#import time
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
	parser = argparse.ArgumentParser()
	#1:integration site, 2:whole reads
	parser.add_argument('-s', dest='separate', help="1: integration site; 2: whole fragment",type=int, default=1)
	parser.add_argument("-shift",help="shift size from integration site(bp), default: 10",dest='shift',default=10,type=int)
	parser.add_argument("-ES",help="extend size from integration site (bp), default: 20",dest='extend',default=20,type=int)
	parser.add_argument("-bs",help="bin size for bigwig (bp), default: 10",dest='binsize',default=10,type=int)
	parser.add_argument('input_bam')
	parser.add_argument('output_name')

	args = parser.parse_args()

	input_bam=args.input_bam
	Outname=args.output_name

	shiftSize=args.shift
	#extend value from intesite
	extSize=args.extend
	bs=args.binsize
	#clean bam files
	#index :print error message (try)

	subprocess.call('''samtools index %s'''%(input_bam),shell=True)

	#option1 integration site/ tn5 cutting site
	if (args.separate == 1):
		# shift: move read 10 ,extend left and right 20
		subprocess.call('''macs2 callpeak -t %s --nomodel --shift -%s --extsize %s -n %s 2>&1'''%(input_bam,shiftSize,extSize,Outname+'.integ_peak'), shell=True)
		#Making bigwig/ bamcoverage limit: cannot use minus
                #subprocess.call('''bamCoverage -b %s -bs 10 --normalizeUsingRPKM --Offset 1 20 -o %s'''%(input_bam, input_bam+'_coverage.bw'),shell=True
                subprocess.call('''bamCoverage -b %s -bs %s --normalizeUsing RPKM --Offset 1 20 -o %s 2>&1'''%(input_bam, bs, Outname+'_coverage.bw'),shell=True)
        elif(args.separate == 2):
        
            subprocess.call('''samtools index %s'''%(input_bam),shell=True)
            #subprocess.call('''samtools index %s'''%(input_bam+'_short.bam'),shell=True)
       
            #-f: paired end
            subprocess.call('''macs2 callpeak -t %s --format BAMPE -n %s 2>&1'''%(input_bam,Outname+'.narrowpeak'), shell=True)
            #subprocess.call('''macs2 callpeak -t %s --format BAMPE --broad -n %s'''%(input_bam+'_short.bam',input_bam+'_short_peak'), shell=True)
            #Making bam coverage
            #subprocess.call('''bamCoverage -b %s -bs 10 --normalizeUsingRPKM --Offset 1 20 -o %s'''%(input_bam, input_bam+'_coverage.bw'),shell=True)
            subprocess.call('''bamCoverage -b %s -bs 10 --normalizeUsing RPKM -e -o %s 2>&1'''%(input_bam, Outname+'_coverage.bw'),shell=True)
            #subprocess.call('''bamCoverage -b %s -bs 10 --normalizeUsingRPKM -e -o %s'''%(input_bam+'_short.bam', input_bam+'_short_coverage.bw'),shell=True)
	
	else:

            print "Please enter separation number"




if __name__ == '__main__':
	main()

