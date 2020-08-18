# ATACgraph

ATACgraph is a simple and effective software for the analysis of ATAC-Seq data. It contains 11 analyses in 9 major modules to profile (epi)genome. 

![ATACgraph flow](https://github.com/RitataLU/ATACgraph_v2/blob/master/ATACgraph%20flow.png)



# <a name="SystemRequirements"></a>System Requirements
* Linux or Mac OS Environment
* Python 2.7
* deepTools 3.3.0
* [BEDtools](http://bedtools.readthedocs.org/) 
* [IDR](https://github.com/nboley/idr)
```
IDR require python 3.5
```

* [MACS2](https://github.com/taoliu/MACS)
* Python Modules:
  * Numpy
  * pandas
  * pysam==0.11.2.2
  * matplotlib
  * matplotlib-venn
  * argparse
  * pybedtools
  * deepTools==3.3.0
  * scikit-learn
  

# Installation
## Docker 
Link: https://hub.docker.com/r/lsbnb/galaxy_atacgraph

## or 
## Linux Command 
1. Download the source code and install the requirements.


```
$ git clone https://github.com/RitataLU/ATACgraph.git
$ cd ATACgraph
$ sudo sh ./base.txt

``` 

2. Add your ATACgraph path to the PATH.

(1)  Edit bash profile
  
``` 
$ vi ~/.bash_profile
``` 
   
   (2) Add ATAC-graph/script path to the PATH environment variable.
 
``` 
$ PATH=$PATH:(ATACgraph/script file path)
$ source ~/.bash_profile

```    

# Running ATACgraph

## ATACgraph modules 


```
$ ATACgraph -h
Usage: atacG <subcommand> [options]
ATACgraph sub-commands include:

00_rmChr                  Remove chrM,chrPt
01_calFragDist            Generate figures of fragments size distribution & Fast Fourier transform
02_selectFragSize         Select bam fragments size
02_gtftoBed               Transform gtf file to bed files
03_callPeak               Call ATAC-seq peaks"
03_genePlot               Generate figures of depicting genes and peaks accessibility
03_junctionBed            Generate junction bed track
04_specificPeaks          Identify specific peaks between 2 groups of peaks
04_specificPeaksIDR       Identify specific peaks with Irreproducibility Discovery Rate (IDR) framework
04_specificPromoter       Identify specific promoter using Gaussian Mixture Model between 2 groups of peaks
05_seqCompare             Compare peaks between ATAC-seq and Other seq
```

## Filtering ATAC-seq reads from any chromosome 

**Input:**
* ATAC-seq reads file in bam

```
$ ATACgraph 00_rmChr -h
Usage: python rmChr.py <input.bam> <output.bam> <chrM>

If you need to remove multiple chromosomes, use comma
to seperate. For example chrM,chrPt

```

## Fragment length distribution and Fast Fourier Transform (FFT)

**Input:**
* ATAC-seq bam file

```
$ ATACgraph 01_calFragDist -h
usage: calFragDist [-h] input_bam output_name

positional arguments:
  input_bam
  output_name

optional arguments:
  -h, --help   show this help message and exit
                        
```

**Output:** 
* 2 figures (fragment lenth distribution & FFT analysis result)

   * fragment distribution
   
   ![fragDist](https://github.com/RitataLU/ATACgraph/blob/master/fragDist.png)
   
   * fragment distribution FFT
   
   ![FFT](https://github.com/RitataLU/ATACgraph/blob/master/FFT.png)


**Output:** 
* filtered bam file 

## Selectiing fragments size 

**Input:**
* ATAC-seq reads file in bam

```
$ ATACgraph 02_selectFragSize -h
usage: 02_selectFragSize [-h] [-f FILTER] [-m MODE] input_bam output_bam

positional arguments:
  input_bam
  output_bam

optional arguments:
  -h, --help            show this help message and exit
  -f FILTER, --filter FILTER
                        default=150
  -m MODE, --mode MODE  Select fragments smaller [1] or larger [2] than filter
                        size. Default=1                       

```

**Output:** 
* fragment bam file 

## Transform GTF file to BED files

**Input:**
* Annotation GTF file

```
$ ATACgraph 02_gtftoBed -h
usage: gtftoBed [-h] [-p PROMOTER] input_gtf gtf_name

positional arguments:
  input_gtf
  gtf_name

optional arguments:
  -h, --help            show this help message and exit
  -p PROMOTER,          --promoter PROMOTER, promoter region from gene transcript start site (TSS) default = 2000
                        
```

**Output:** 
* 8 BED file (promoter,gene,exon,intron,utr5,cds,utr3,igr) for the generating metagene plots, fold enrichment analysis 



## Generating fragment size tracks 

**Input:**
* ATAC-seq bam file

```
ATACgraph 03_junctionBed -h
usage: junctionBed [-h] [-s SEPARATE] [-b BIN] [-f FILTER]
                   input_bam output_bed

positional arguments:
  input_bam
  output_bed

optional arguments:
  -h, --help            show this help message and exit:w
  -s SEPARATE, --separate SEPARATE
                        border lenth of long and short fragment(bp), default: 150
  -b BIN, --bin BIN     binzise of genome(bp), default: 10
  -f FILTER, --filter FILTER
                        number of fragment juction tracks, default:1

```

**Output:** 

* track BED files
* Visualization on IGV

![specificPeak](https://github.com/RitataLU/ATACgraph/blob/master/Figures/Fig3A.png)

##  ATAC-seq peak calling

**Input:**
* ATAC-seq bam file

```
ATACgraph 03_callPeak -h
usage: callPeak [-h] [-s SEPARATE] [-shift SHIFT] [-ES EXTEND] [-bs BINSIZE]
                input_bam output_name gene_bed

positional arguments:
  input_bam     input bam file
  output_name   name for output files
  gene_bed      Gene or promoter annotation bed file, either
                gene_body_bed6.bed or gene_promoter_bed6.bed

optional arguments:
  -h, --help    show this help message and exit
  -s SEPARATE   1: integration site; 2: full-extend fragment
  -shift SHIFT  shift size from integration site(bp), default: 50
  -ES EXTEND    extend size from integration site (bp), default: 100
  -bs BINSIZE   bin size for bigwig (bp), default: 10
  
```

**Output:** 
* Peak location BED file (.narrowpeak), 
* Peak intensity bigWigfile (.coverage.bw) 
* A genes list of overlapping with peaks locations


## Identify differential enriched ATAC-seq peaks between two conditions.

To identify differentially accessible regions, ATACgraph utilizes peak BED files and BigWig files generated by callPeak module to compute peaks abundances between two samples.

**Input:**
* 2 ATAC-seq peaks BDE files
* 2 ATAC-seq abundandance BigWig files


```
ATACgraph 04_specificPeaks -h
usage: specificPeaks [-h] [-c FOLD_CHANGE] [-p P_VALUE] input_peakAs input_peakBs input_bwAs input_bwBs output_bedA output_bedB

positional arguments:
  input_peakAs  Peak beds, seprate by comma
  input_peakBs  Peak beds, seprate by comma
  input_bwAs    BigWig files, seprate by comma
  input_bwBs    BigWig files, seprate by comma
  output_bedA
  output_bedB

optional arguments:
  -h, --help    show this help message and exit
  -c FOLD_CHANGE, --fold_change FOLD_CHANGE
                        Fold change cutoff. Default:2
  -p P_VALUE, --p_value P_VALUE
                        P-value cutoff. Default:0.05
  
```
**Output:** 

* 2 ATAC-seq specific peaks BED files


## Identify differential enriched ATAC-seq peaks between two conditions with IDR

```
idr --samples 1. narrowPeak 2.narrowPeak --output-file idr.txt
```
**Output:** 

* A txt files, format information are based on [IDR](https://github.com/nboley/idr)

```
ATACgraph 04_specificPeaks -h
usage: specificPeaks [-h] [-c FOLD_CHANGE] [-p P_VALUE] input_peakAs input_peakBs output_bedA output_bedB

positional arguments:
  input_peakAs  Peak beds, seprate by comma
  input_peakBs  Peak beds, seprate by comma
  output_bedA
  output_bedB

optional arguments:
  -h, --help    show this help message and exit
  
```
**Output:** 

* 2 ATAC-seq specific peaks BED files


## Comparing peaks between ATAC-seq and another Seq.

**Input:**
* ATAC-seq peaks BED file & another seq peaks file

```
ATACgraph  05_seqCompare -h
usage: seqCompare [-h] atac_peak other_peak ATAC_name otherPeak_name overlap_name Genome_size genes

positional arguments:
  atac_peak       ATAC-seq peak bed
  other_peak      Other seq peak bed
  ATAC_name       Name for ATAC peak
  otherPeak_name  Name for other peak
  overlap_name    Name for overlapping peak
  Genome_size     Genome size(bp)
  genes           Gene or promoter annotation file with bed6.bed format

optional arguments:
  -h, --help      show this help message and exit
```
**Output:** 
* An overlapping peaks location BED file 
* Venn diagram shows the numbers of each seq peaks and over lapping peaks with p value (hypergeometric test)
* A gene list of overlapping peaks locations between 2 seq

![Venn](https://github.com/RitataLU/ATACgraph/blob/master/venn.png)


## Heatmap and metagene plots of ATAC-seq abundance, Fold enrichment analysis of open regions in genomic features 
To investigate the chromatin accessibility around genes, To investigate the chromatin accessibility around genes, ATACgraph uses the files describing the ATAC-seq peak locations and gene annotations for two types of analyse. This step requires 8 genomic feature BED files, user should run gtftoBed before this step.

**Input:**
* ATAC-seq bam file

```
ATACgraph 03_genePlot -h
usage: genePlot [-h] [-p PROMOTER] input_peak input_bigwig gtf_name

positional arguments:
  input_peak
  input_bigwig
  gtf_name

optional arguments:
  -h, --help            show this help message and exit
  -p PROMOTER, --promoter PROMOTER
  

```

**Output:** 
*  text files
   * value of Heatmap depicting accessibility for gene (genebody.matrix.txt & genebody.matrix.gz)
   * value of Heatmap depicting accessibility for peak (peak.matrix.txt & peak.matrix.gz)
   * The intersection site between 8 genomic features and peaks (8 files)
   
* Figures 
   * Fold enrichment analysis of open regions in genomic features (.Fold_Enrichment.png)
   ![Enrichment](https://github.com/RitataLU/ATACgraph_v2/blob/master/FoldEnrichment.png)
   
   
   * heatmap, metaplot ATAC-seq abundance related to genes     
   * heatmap, metaplot ATAC-seq abundance related toand peaks
  
   <img align="left" width="300" height="900" src="https://github.com/RitataLU/ATACgraph/blob/master/TKO.integ_coverage.bwgene_body_heatmap.png">

   <img align="right" width="300" height="900" src="https://github.com/RitataLU/ATACgraph/blob/master/Peak_heatmap.png">
   



