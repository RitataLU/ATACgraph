
#20201204 add a parameter for control bam file
# callpeak & make bw file seperate only integration, no integration
import sys
import pandas as pd 
import os
import subprocess, sys
#import time
import argparse

def main():
	parser = argparse.ArgumentParser()
	#1:integration site, 2:whole reads
	parser.add_argument('-s', dest='separate', help="1: integration site; 2: full-extend fragment",type=int, default=2)
	parser.add_argument("-shift",help="shift size from integration site(bp), default: 50",dest='shift',default=50,type=int)
	parser.add_argument("-ES",help="extend size from integration site (bp), default: 100",dest='extend',default=100,type=int)
	parser.add_argument("-bs",help="bin size for bigwig (bp), default: 10",dest='binsize',default=10,type=int)
	parser.add_argument("-c",dest = "control_bam",help='input control bam file, default: none',default='',type=str)
	parser.add_argument('input_bam',help='input ATAC-seq bam file',type=str)
	parser.add_argument('output_name',help='name for output files',type=str)
	parser.add_argument('gene_bed',help='Gene or promoter annotation bed file, either gene_body_bed6.bed or gene_promoter_bed6.bed',type=str)

	args = parser.parse_args()

	input_bam=args.input_bam
	Outname=args.output_name
	control_bam = args.control_bam 
	control_para = '' if not control_bam else '-c '+control_bam
	shiftSize=args.shift
	#extend value from intesite
	extSize=args.extend
	bs=args.binsize
	gene = args.gene_bed
        #clean bam files
	#index :print error message (try)

	subprocess.call('''samtools index %s'''%(input_bam),shell=True)

	#option1 integration site/ tn5 cutting site
	if (args.separate == 1):
		# shift: move read -50 ,extend left and right 100
		subprocess.call('''macs2 callpeak %s -t %s --nomodel --shift -%s --extsize %s -n %s 2>&1'''%(control_para,input_bam,shiftSize,extSize,Outname), shell=True)
		#Making bigwig/ bamcoverage limit: cannot use minus
		#subprocess.call('''bamCoverage -b %s -bs %s --normalizeUsingRPKM --Offset 1 20 -o %s 2>&1'''%(input_bam, input_bam+'_coverage.bw'),shell=True
		subprocess.call('''bamCoverage -b %s -bs %s --normalizeUsing RPKM --Offset 1 20 -o %s 2>&1'''%(input_bam, bs, Outname+'_coverage.bw'),shell=True)
		#subprocess.call('''bedtools -a %s -b %s -wa -wo '''%(Outname+'.peak',gene,Outname+'peak_Gene_list.txt'),shell=True)
		subprocess.call('''bedtools intersect -wo -a %s -b %s >%s '''%(gene,Outname+'_peaks.narrowPeak',Outname+'_peak_gene_list.txt'),shell=True)

	elif (args.separate == 2):
		#-f: paired end
		subprocess.call('''macs2 callpeak %s -t %s --format BAMPE -n %s 2>&1'''%(control_para,input_bam,Outname), shell=True)
		#Making bam coverage
		#subprocess.call('''bamCoverage -b %s -bs %s --normalizeUsing RPKM -e -o %s 2>&1'''%(input_bam, bs, Outname+'_coverage.bw'),shell=True)
		subprocess.call('''bamCoverage -b %s -bs %s --normalizeUsing RPKM -e -o %s 2>&1'''%(input_bam, bs, Outname+'_coverage.bw'),shell=True)
		#subprocess.call('''bedtools -a %s -b %s -wa -wo '''%(Outname+'.peak',gene,Outname+'peak_Gene_list.txt'),shell=True)
		subprocess.call('''bedtools intersect -wo -a %s -b %s >%s '''%(gene,Outname+'_peaks.narrowPeak',Outname+'_peak_gene_list.txt'),shell=True)

	else:

		print ("Please enter separation number, 1:")

		subprocess.call('''bedtools intersect -wo -a %s -b %s >%s '''%(gene,Outname+'_peaks.narrowPeak',Outname+'_peak_gene_list.txt'),shell=True)

        
		g1 = pd.read_csv(Outname+"_peak_gene_list.txt",sep="\t",header=None)
		g1.columns=['chr','sta','end','gene','score','dir','chrp','sta_p','end_p','peakname','score_p','strand', 'signalValue','pval','qValue','peak','overlap']


		g1['Peak']="("+g1['chrp']+":"+g1['sta_p'].astype(str)+"-"+g1['end_p'].astype(str)+')'

		g2 = g1[['chr','sta','end','gene','Peak']]
		g2.to_csv(Outname+'_peak_gene_list.txt',sep="\t",index = None)




if __name__ == '__main__':
	main()

