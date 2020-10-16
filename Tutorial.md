# Tutorial

# Example use case 

## Installation
### Download the source code and install the requirements.

```
$ git clone https://github.com/RitataLU/ATACgraph.git
$ cd ATACgraph
$ sudo sh ./base.txt

``` 

### Add your ATACgraph path to the PATH.

(1)  Edit bash profile
  
``` 
$ vi ~/.bash_profile
``` 
   
(2) Add ATAC-graph/script path to the PATH environment variable.
 
``` 
$ PATH=$PATH:(ATACgraph/script file path)
$ source ~/.bash_profile

```    
## Run demo
Download the demo input file in ATACgraph folder

```
$ cd ATACgraph/demo
$ tar -xvf demo.tar.gz 
``` 
### Remove mitochondria chrmosome

**Input:**
* ATAC-seq bam file
``` 
$ ATACgraph 00_rmChr demo.bam demo_rmM.bam chrM
```
**Output:**
* ATAC-seq bam file after removing mitochondria chromosome named demo_rmM.bam

## Fragment length distribution and Fast Fourier Transform (FFT)

**Input:**
* ATAC-seq bam file after removing mitochondria chromosome

```
$ ATACgraph 01_calFragDist demo_rmM.bam demo_rmM_fragment demo_rmM_FFT

```
              
**Output:** 
* 2 figures (demo_rmM_fragment.png & demo_rmM_FFT.png)


## Selectiing fragments size 

**Input:**
* ATAC-seq bam file after removing mitochondria chromosome

```
$ samtools index demo_rmM.bam
$ ATACgraph 02_selectFragSize demo_rmM.bam demo_rmM_output.bam [-f 150] [-m 1]

optional arguments:
  -h, --help            show this help message and exit
  -f FILTER, --filter FILTER
                        default=150
  -m MODE, --mode MODE  Select fragments smaller [1] or larger [2] than filter
                        size. Default=1                       

```
**Output:** 
* demo_rmM_output.bam


## Transform GTF file to BED files

**Input:**
* Annotation GTF file

```
$ ATACgraph 02_gtftoBed demo_gene.gtf demo [-p 2000]

optional arguments:
  -h, --help            show this help message and exit
  -p PROMOTER,          --promoter PROMOTER, promoter region from gene transcript start site (TSS) default = 2000
                        
```

**Output:** 
* 11 BED file (promoter,gene,exon,intron,utr5,cds,utr3,igr) for the generating metagene plots, fold enrichment analysis 



## Generating fragment size tracks 

**Input:**
* ATAC-seq bam file after removing mitochondria chromosome

```
ATACgraph 03_junctionBed demo_rmM.bam demo_rmM_junction_output.bed

```

**Output:** 

* A track BED file named demo_rmM_junction_output.bed for visualizing on IGV


##  ATAC-seq peak calling

**Input:**
* ATAC-seq bam file after removing mitochondria chromosome

```
ATACgraph 03_callPeak demo_rmM.bam demo_rmM_peakcall demo_gene_body_bed6.bed
  
```

**Output:** 
* Peak location BED file (demo_rmM_peakcall_peaks.narrowPeak), 
* Peak intensity bigWigfile (demo_rmM_peakcall_coverage.bw) 
* A genes list of overlapping with peaks locations (demo_rmM_peakcall_peak_gene_list.txt)




## Heatmap and metagene plots of ATAC-seq abundance, Fold enrichment analysis of open regions in genomic features 
To investigate the chromatin accessibility around genes, To investigate the chromatin accessibility around genes, ATACgraph uses the files describing the ATAC-seq peak locations and gene annotations for two types of analyse. This step requires 8 genomic feature BED files, user should run gtftoBed before this step.

**Input:**
* ATAC-seq bam file after removing mitochondria chromosome 

```
ATACgraph 03_genePlot demo_rmM_peakcall.narrowpeak demo_rmM_peakcall_coverage.bw demo 
  
```

**Output:** 
* 3 Figures
  * The enrichment status of accessible region in genome (Fold_Enrichment.pdf)
  * The accessibility – or read abundance – around genes (gene_body_heatmap.pdf)
  * The accessibility – or read abundance – around peaks (Peak_heatmap.pdf)
  
*  text files
   * value of Heatmap depicting accessibility for gene (demo_rmM_peakcall_coverage.bwgene_body.matrix.txt & demo_rmM_peakcall_coverage.bwgene_body.matrix.gz)
   * value of Heatmap depicting accessibility for peak (demo_rmM_peakcall_coverage.bw_peak.matrix.txt & demo_rmM_peakcall_coverage.bw_peak.matrix.gz)
   * The value to generate enrichment figure(demo_rmM_peakcall_peaks.narrowPeak_Fole_Enrichment_Table.txt)  
   * The intersection site between 8 genomic features and peaks: 
     * demo_rmM_peakcall_peaks.narrowPeak_3utr.txt
     * demo_rmM_peakcall_peaks.narrowPeak_exons.txt                  
     * demo_rmM_peakcall_peaks.narrowPeak_gene_igr.txt
     * demo_rmM_peakcall_peaks.narrowPeak_5utr.txt         
     * demo_rmM_peakcall_peaks.narrowPeak_gene_promoter.txt
     * demo_rmM_peakcall_peaks.narrowPeak_cds.txt   
     * demo_rmM_peakcall_peaks.narrowPeak_gene_body.txt             
     * demo_rmM_peakcall_peaks.narrowPeak_introns.txt




