# ATACgraph


* Python 2.7

    `SAMtools <http://www.htslib.org/>`_ 
    `deepTools <https://deeptools.readthedocs.org>`_
    `BEDtools <http://bedtools.readthedocs.org/>`_ 

* Python Modules 'Numpy', 'pandas' and 'Metplotlib'. To install the packages, use the following commands on an UNIX terminal:
 
    $ pip install numpy
    $ pip install pandas
    $ pip install matplolib
  
  Installation
============

1. Download the source code and install the requirements.

  ::

  $ git clone https://github.com/kullatnunu/atacgraph.git
  $ pip install -r atacgraph/base.txt

  
2. Add your ATAC-graph path to the PATH.

   (1) Edit bash profile
  
   ::
  
   $ vi ~/.bash_profile
   
   (2) Add ATAC-graph path to the PATH environment variable.
 
   ::
  
   $ PATH=$PATH:(ATAC-graph file path)
   $ source ~/.bash_profile
   
Running ATAC-graph
==================
