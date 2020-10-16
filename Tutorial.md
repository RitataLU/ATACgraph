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
### remove mitochondria chrmosome

``` 
$ ATACgraph 00_rmChr hESC_TKO_chr22.bam hESC_TKO_chr22_rmM.bam chrM
```

## Fragment length distribution and Fast Fourier Transform (FFT)

