# ATACgraph

ATACgraph is a simple and effective software for the analysis of ATAC-Seq data. It contains 11 analyses in 9 major modules to profile (epi)genome. 

![ATACgraph flow](https://github.com/RitataLU/ATACgraph_v2/blob/master/ATACgraph%20flow.png)



# <a name="SystemRequirements"></a>System Requirements
* Linux or Mac OS Environment
* Python2.7 (it should be pre-installed in both Linux and Mac). Type 'Python' to see the installed version. Python2 could be downloaded from http://www.python.org/download )
* [deepTools](https://deeptools.readthedocs.org/)
* [BEDtools](http://bedtools.readthedocs.org/) 
* Python Modules. To install the packages, use the following commands on an UNIX terminal:


```
pip install numpy
pip install pandas
pip install pysam==0.11.2.2
pip install matplotlib
pip install argparse
pip install pybedtools

```

# Installation


1. Download the source code and install the requirements.


```
$ git clone https://github.com/RitataLU/ATACgraph.git
$ pip install -r ATACgraph/.txt

``` 

2. Add your ATACgraph path to the PATH.

   (1) Edit bash profile
  
``` 
$ vi ~/.bash_profile
``` 
   
   (2) Add ATAC-graph path to the PATH environment variable.
 
``` 
$ PATH=$PATH:(ATAC-graph file path)
$ source ~/.bash_profile
```    

# Running ATACgraph

## ATACgraph modules 

```
$ ATACgraph -h
Usage: atacG <subcommand> [options]
atacG sub-commands include:

rmChr             Remove chrM,chrPt
calFragDist       Generate figures of fragments size distribution & Fast Fourier transform
selectFragSize    Select bam fragments size
junctionBed       Generate junction bed track
gtftoBed          Transform gtf file to bed files
callPeak          call ATAC-seq peaks
specificPeaks     Find specific peaks between 2 peaks files
genePlot          Generate figures of depicting gene accessibility
seqCompare        Compare peaks between ATAC-seq and Other seq
```

## Filtering ATAC-seq reads from any chromosome 

**Input:**
* ATAC-seq reads file in bam

```
$ ATACgraph rmChr -h
Usage: python rmChr.py <input.bam> <output.bam> <chrM>

If you need to remove multiple chromosomes, use comma
to seperate. For example chrM,chrPt

```

**Output:** 
* filtered bam file 

## Fragments size selection 

**Input:**
* ATAC-seq reads file in bam

```
$ ATACgraph selectFragSize -h
usage: selectFragSize [-h] [-f FILTER] [-m MODE] input_bam output_bam

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
$ ATACgraph gtftoBed -h
usage: gtftoBed [-h] [-p PROMOTER] input_gtf gtf_name

positional arguments:
  input_gtf
  gtf_name

optional arguments:
  -h, --help            show this help message and exit
  -p PROMOTER,          --promoter PROMOTER, promoter region from gene transcript start site (TSS) default = 2000
                        
```

**Output:** 
* 8 BED file (promoter,gene,exon,intron,utr5,cds,utr3,igr)


## Fragment length distribution and Fast Fourier Transform (FFT)

**Input:**
* ATAC-seq bam file

```
$ ATACgraph calFragDist -h
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
   
   ![fragDist](https://github.com/RitataLU/ATACgraph_v2/blob/master/fragDist.png)
   
   * fragment distribution FFT
   
   ![FFT](https://github.com/RitataLU/ATACgraph_v2/blob/master/FFT.png)


##  ATAC-seq peak calling

**Input:**
* ATAC-seq bam file

```
ATACgraph callPeak -h
usage: callPeak [-h] [-s SEPARATE] [-shift SHIFT] [-ES EXTEND] [-bs BINSIZE]
                input_bam output_name

positional arguments:
  input_bam
  output_name

optional arguments:
  -h, --help    show this help message and exit
  -s SEPARATE   1: integration site; 2: whole fragment
  -shift SHIFT  shift from integration site, default=10
  -ES EXTEND    extend size from integration site, default=20
  -bs BINSIZE   binsize for bigwig, default=10
  
```

**Output:** 
* Peak location BED file (.narrowpeak) & Peak intensity bigWigfile (.coverage.bw)


## Heatmap and metagene plots of ATAC-seq abundance, Fold enrichment analysis of open regions in genomic features 

**Input:**
* ATAC-seq bam file

```
ATACgraph genePlot -h
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

* Figures 
   * heatmap, metaplot ATAC-seq abundance related to genes 
   * heatmap, metaplot ATAC-seq abundance related toand peaks
   
   ![heatmap](https://github.com/RitataLU/ATACgraph_v2/blob/master/TKO.integ_coverage.bwgene_body_heatmap.png)
   
   * Fold enrichment analysis of open regions in genomic features (.Fold_Enrichment.png)
   
   ![Enrichment](https://github.com/RitataLU/ATACgraph_v2/blob/master/FoldEnrichment.png)

*  text files
   * value of Heatmap depicting accessibility for gene (genebody.matrix.txt & genebody.matrix.gz)
   * value of Heatmap depicting accessibility for peak (peak.matrix.txt & peak.matrix.gz)
   * The intersection site between 8 genomic features and peaks (8 files)
   
  
## Generate fragment size tracks 

**Input:**
* ATAC-seq bam file

```
ATACgraph junctionBed -h
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

![specificPeak](https://github.com/RitataLU/ATACgraph_v2/blob/master/junction.png)


## Identify differential enriched ATAC-seq peaks between two conditions.

**Input:**
* 2 ATAC-seq peaks

```
ATACgraph specificPeaks -h
usage: specificPeaks [-h] input_peakAs input_peakBs output_bedA output_bedB

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
ATACgraph  seqCompare -h
usage: seqCompare [-h] atac_peak other_peak ATAC_name otherPeak_name overPeak

positional arguments:
  atac_peak       ATAC-seq peak bed
  other_peak      Other seq peak bed
  ATAC_name       name for ATAC peak
  otherPeak_name  name for other peak
  overPeak        name for overlapping peak

optional arguments:
  -h, --help      show this help message and exit
```
**Output:** 

* Venn diagram

![Venn](https://github.com/RitataLU/ATACgraph_v2/blob/master/venn.png)


