# ATACgraph

ATACgraph is a simple and effective software for the analysis of ATAC-Seq data. It contains 11 analyses in 9 major modules to profile (epi)genome. 

# <a name="SystemRequirements"></a>System Requirements
* Linux or Mac OS Environment
* Python2.7 (it should be pre-installed in both Linux and Mac). Type 'Python' to see the installed version. Python2 could be downloaded from http://www.python.org/download/ )

* Python Modules 'Pysam' and 'Metplotlib'. To install the packages, use the following commands on an UNIX terminal:

```
pip install numpy
pip install pandas
pip install pysam
pip install Matplolib

```

# Installation


1. Download the source code and install the requirements.


```
$ git clone https://github.com/kullatnunu/atacgraph.git
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

# Tutorial

##1. remove mithochdria : rmChr 
