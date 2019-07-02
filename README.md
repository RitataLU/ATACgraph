# ATACgraph

ATACgraph is a simple and effective software for the analysis of ATAC-Seq data. It contains 11 analyses in 9 major modules to profile (epi)genome. 

![ATACgraph flow](https://github.com/RitataLU/ATACgraph_v2/blob/master/ATACgraph%20flow.png)



# <a name="SystemRequirements"></a>System Requirements
* Linux or Mac OS Environment
* Python2.7 (it should be pre-installed in both Linux and Mac). Type 'Python' to see the installed version. Python2 could be downloaded from http://www.python.org/download )

* To install the packages, use the following commands on an UNIX terminal:

```
pip install numpy
pip install pandas
pip install pysam
pip install Matplolib
pip install argparse
pip install pybedtools

```

# Installation


1. Download the source code and install the requirements.


```
$ git clone https://github.com/RitataLU/atacgraph_v2.git
$ pip install -r atacgraph/base.txt

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
$ TACgraph selectFragSize -h
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
  -p PROMOTER, --promoter PROMOTER
                        
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
   * fragment_distribution.png
   * fragment_distribution_fft.png










