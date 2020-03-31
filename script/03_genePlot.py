
#making , fold enrichment, heatmap
import pandas as pd 
import sys
import os
import numpy as np
import subprocess, sys
import time
import math
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl

#output: heatmap.png, enrichment.png, matrix.gz, .matrix.txt', peak & genomic overlap.txt
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
                

def genome_size(samplename):
        genome_bp=pd.read_csv(samplename+'_genome.bed', sep='\t', header=None)
        genome_bp = genome_bp[1].sum()
        genome_bp=genome_bp*1.0
        return genome_bp


#peak size count
def peak_bp(peak):
        peak = pd.read_csv(peak, header=None, sep="\t").ix[:,0:5]
        peak.columns = ['chr', 'peak_str', 'peak_end', 'peak_id', 'peak_value', 'peak_dir']
        peak=peak.drop_duplicates(subset=['peak_str', 'peak_end'], keep='first')
        peak['bp'] = (peak.peak_end - peak.peak_str)
        peak_bp = peak.bp.sum()
        return peak_bp


#annotation size count
def ano_bp(anno,samplename):
	#anno=genebody, utr...
        ano = pd.read_csv(samplename+'_'+anno+'_merge.bed', header=None, sep="\t")
        ano.columns=['chr','g_str','g_end','gene_id','g_id','g_dir']
        ano=ano.drop_duplicates(subset=['g_str', 'g_end'], keep='first')
        ano['bp'] = (ano.g_end - ano.g_str)
        ano_bp = ano.bp.sum()
        return ano_bp

#annotation_peak bp count
def annopeakbp(peak,anno):
        anno_peak_asso = pd.read_csv(peak+'_'+anno+'.txt', header=None, sep="\t")
        anno_peak_asso_bp = anno_peak_asso[len(anno_peak_asso.columns)-1].sum()
        anno_peak_asso_bp = float(anno_peak_asso_bp)
        return anno_peak_asso_bp



#making Enrichment Graph
def enrichment_num(peak,anno,genome_bp,samplename):
        try:
                enrichment=math.log((annopeakbp(peak,anno)/peak_bp(peak))/(ano_bp(anno,samplename)/genome_bp),2)
                return enrichment
        except pd.errors.EmptyDataError:
                return 0


#making associate table
#making table and peak_bed to bedgraph

def associate(input_peak,samplename,promoter):
        #hg19testMain_gene_body_merge.bed
        peak = pd.read_csv(input_peak, header=None, sep="\t").ix[:,0:5]
        gene_body = pd.read_csv(samplename+'_gene_body_merge.bed', header=None, sep="\t")
        peak.columns = ['chr', 'peak_str', 'peak_end', 'peak_name','peak_value','peak_dir']
        peak.chr=peak.chr.astype(str)
        gene_body.columns = ['chr', 'gbed_str', 'gbed_end', 'gene_id', 'gene_value', "gene_dir"]
        gene_body['pro_str'] = np.where(gene_body.gene_dir == '+', gene_body.gbed_str - promoter, gene_body.gbed_end -0)
        gene_body['pro_end'] = np.where(gene_body.gene_dir == '+', gene_body.gbed_str +0, gene_body.gbed_end + promoter)
        combined = pd.merge(gene_body, peak, on='chr')
        combined['Genebody'] = np.where((combined.gbed_str < combined.peak_end) & (combined.gbed_end > combined.peak_str), 1, 0)
        combined['Promoter'] = np.where((combined.pro_str < combined.peak_end) & (combined.pro_end > combined.peak_str), 1, 0)
        summary = combined[(combined.Promoter > 0) | (combined.Genebody > 0)]
        s1 = summary.drop(summary.columns[6:13], axis = 1).drop(summary.columns[4], axis = 1)
        s1_group=s1.groupby(['gene_id', 'chr', 'gbed_str', 'gbed_end', 'gene_dir']).agg({'Genebody':'sum', 'Promoter':'sum'}).reset_index().sort_values(["chr","gbed_str"])
        s1_group.to_csv(input_peak+"_gene"+"_summary_table.xls", sep='\t',index=False)

def coverage_heatmap(coverage,samplename,input_peak):
        subprocess.call('''computeMatrix scale-regions -S %s -R %s --missingDataAsZero -bs 10 -a 1000 -b 1000 -out %s --outFileNameMatrix %s 2>&1'''%(coverage,samplename+'_gene_body_merge.bed',coverage+'gene_body'+'.matrix.gz',coverage+'gene_body'+'.matrix.txt'),shell=True)
        subprocess.call('''plotHeatmap -m %s -out %s --legendLocation none'''%(coverage+'gene_body'+'.matrix.gz','gene_body_heatmap.pdf'),shell=True)

        #peak heatmap
        subprocess.call('''computeMatrix reference-point --referencePoint center -S %s -R %s --missingDataAsZero -bs 10 -a 1000 -b 1000 -out %s --outFileNameMatrix %s 2>&1'''%(coverage,input_peak,coverage+'_peak'+'.matrix.gz',coverage+'_peak'+'.matrix.txt'),shell=True)

        subprocess.call('''plotHeatmap --refPointLabel peak  --regionsLabel peaks --xAxisLabel 'peak distance(bp)' -m %s -out %s --legendLocation none'''%(coverage+'_peak'+'.matrix.gz','Peak_heatmap.pdf'),shell=True)

#Making annotation_peak_associate file
def annotation_peak(input_peak,samplename):
        gene = "gene_body"
        exon = "exons"
        intron = "introns"
        utr3 = "3utr"
        utr5 = "5utr"
        cds = "cds"
        promoter = "gene_promoter"
        igr = "gene_igr"
        annotation_name=[promoter,gene,exon,intron,utr5,cds,utr3,igr]

        file_exist=[]
        for i in annotation_name:
                subprocess.call('''bedtools intersect -nonamecheck -a %s -b %s -wo > %s'''%(samplename+'_'+i+'_merge.bed',input_peak,input_peak+'_'+i+'.txt'),shell=True)
        #subprocess.call('''bedtools intersect -nonamecheck -a %s -b %s -wo > %s'''%(samplename+'_gene_body_merge.bed',input_peak,input_peak+'_gene_body.txt'),shell=True)
                if os.stat( input_peak+'_'+i+'.txt' ).st_size == 0 :
                        file_exist.append(1)

def enrichment(input_peak,genome_bp,samplename):

        gene = "gene_body"
        exon = "exons"
        intron = "introns"
        utr3 = "3utr"
        utr5 = "5utr"
        cds = "cds"
        promoter = "gene_promoter"
        igr = "gene_igr"
        annotation_name=[promoter,gene,exon,intron,utr5,cds,utr3,igr]      
 
        #making fold enrichment graph

        plt.style.use('ggplot')
        annotationname = ['Promoter','Genebody','Exon','Intron','5UTR','CDS','3UTR','IGR']
        annotationname_index = range(len(annotationname))
        enrichment_data = []
        for j in annotation_name:
                enrichment_data.append(enrichment_num(input_peak,j,genome_bp,samplename))

        fold_enrich = enrichment_data
        fig = plt.figure()
        ax1 = fig.add_subplot(1,1,1)
        colors=['firebrick','darkblue','gold','darkviolet','darkgreen','darkcyan','c','purple']
        ax1.bar(annotationname_index, fold_enrich, align='center',color=colors)
        ax1.xaxis.set_ticks_position('bottom')
        ax1.yaxis.set_ticks_position('left')
        plt.ylim(ymin=-2,ymax=6)
        plt.xticks(annotationname_index, annotationname, fontsize='small')
        plt.ylabel("Fold Enrichment (log2)")
        plt.title('')
        plt.savefig('Fold_Enrichment.pdf',dpi=400,bbox_inches='tight')
        plt.close(fig)
        fe_table=pd.DataFrame([fold_enrich],columns=annotationname)
        fe_table.to_csv(input_peak+'_Fole_Enrichment_Table.txt', index=None, sep="\t")

#def main():
    #if len(sys.argv) != 4:
        #print "Usage: python rmChr.py <input.bam> <output.bam> <chrM>"
        #print ""
        #print "If you need to remove multiple chromosomes, use comma"
        #print "to seperate. For example chrM,chrPt"
        #sys.exit()

#    infile = 'Ctrl_1.bam'
#    outfile = "test.bam"
#    targetChr = 'chr1,chrX,chrY,chr20'
    #infile = sys.argv[1]
    #outfile = sys.argv[2]
    #targetChr = sys.argv[3]
    #rmChr(infile,outfile,targetChr)

def main():

        parser = argparse.ArgumentParser()
        parser.add_argument('input_peak', type=str)
        parser.add_argument('input_bigwig', type=str)
        parser.add_argument('gtf_name', type=str)
        parser.add_argument('-p', '--promoter', type=int, default=2000)
        args = parser.parse_args()

        promoter=args.promoter
        input_peak=args.input_peak
        samplename=args.gtf_name
        bigwig=args.input_bigwig

        summary_table = associate(input_peak,samplename,promoter)
        bamcoveragegraph=coverage_heatmap(bigwig,samplename,input_peak)
        annotation_peak(input_peak,samplename)
        genome_bp=genome_size(samplename)
        enrichment(input_peak,genome_bp,samplename)



#subprocess.call('''rm *bed6*|rm *3utr.bed|rm *5utr.bed|rm *_start.bed|rm *_cds.bed|rm *_codingexons.bed|rm *_codingintrons.bed|rm *_exons.bed|rm *_igr.bed|rm *_genome.bed|rm *_introns.bed|rm *_noncodingexons.bed|rm *_noncodingintrons.bed|rm *bam.bai|rm *sam|rm *bam_chr|rm *gappedPeak|rm *peaks.xls|rm *.bampaired ''',shell=True)

if __name__ == '__main__':
        main()



 #"***----------Processing Time: %s seconds ----------***" %(tend-tstart)

