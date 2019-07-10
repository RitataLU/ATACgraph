#compare with DNase-seq : postive control
import pandas as pd
import subprocess, sys
import time
import argparse
import subprocess
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

def main():
	parser = argparse.ArgumentParser()
        parser.add_argument('atac_peak',help='ATAC-seq peak bed',type=str)
	parser.add_argument('other_peak',help='Other seq peak bed',type=str)
	parser.add_argument('ATAC_name',help='name for ATAC peak',type=str)
	parser.add_argument('otherPeak_name',help='name for other peak',type=str)
	parser.add_argument('overPeak',help='name for overlap peak',type=str)
	args = parser.parse_args()

	A=args.atac_peak
	D=args.other_peak
	#propotion=args.propotion
	overPeak=args.overPeak
	
	#overlaping bed
	#subprocess.call('''bedtools intersect -a %s -b %s -f %s -wo > %s'''%(A,D,propotion,overPeak+'_'+str(propotion)+'.bed'), shell=True)
	subprocess.call('''bedtools intersect -a %s -b %s -wo > %s'''%(A,D,overPeak+'.bed'), shell=True)
	
	o=str(overPeak)+'.bed'
	atac_peak=pd.read_csv(A,sep="\t",header=None)
	DNase_peak=pd.read_csv(D,sep="\t",header=None)
	oPeak=pd.read_csv(o,sep="\t",header=None)

	# of atac peaks
	a=len(atac_peak)
	# of DNase peaks
	b=len(DNase_peak)
	# overlap peaks
	c=len(oPeak)
	
	otherPeak_name=args.otherPeak_name
	ATAC_name=args.ATAC_name

	out = venn2(subsets = (a-c, b-c, c), set_labels = (ATAC_name,otherPeak_name))
	#out = venn2(subsets = (1,2 , 3), set_labels = ("A","b"))
	out.get_patch_by_id('10').set_color('pink')
	out.get_patch_by_id('01').set_color('peru')
	out.get_patch_by_id('10').set_edgecolor('none')
	out.get_patch_by_id('01').set_edgecolor('none')
	out.get_patch_by_id('11').set_color('olive')
	out.get_patch_by_id('11').set_edgecolor('none')

	for text in out.set_labels:
		text.set_fontsize(14)
	for text in out.subset_labels:
        	text.set_fontsize(14)

	#plt.show()
	plt.savefig(overPeak+'.png')
        

if __name__ == '__main__':
    main()






