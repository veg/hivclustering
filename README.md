HIVClustering
-------------

A Python 3 library that makes infers molecular transmission networks from sequence data. A part of [HIV-TRACE](https://academic.oup.com/mbe/article/35/7/1812/4833215), available at http://hivtrace.hyphy.org

Related projects

1. [BioExt](https://github.com/veg/bioext) Provides the `bealign` utility for codon-aware mapping to a reference sequence.
2. [TN93](https://github.com/veg/tn93) Provides the `tn93` tool for rapid computation of pairwise genetic distances between aligned sequence. 
 
Dependencies
------------

HIVClustering can be installed in one easy step.  Prior to installation, please
ensure that you `python 3` installed. 

`pip3 install hivclustering`

To install a version that has not yet been published, i.e. the `dev` branch, use 

```
git clone https://github.com/veg/hivclustering.git
git checkout develop
python3 setup.py install
```
    
Components
-----

1. [hivnetworkcsv](https://github.com/veg/hivclustering/wiki/hivnetworkcsv): the workhorse of the package which consumes pairwise distance files (and, optionally, other sources of data) and infers the molecular transmission network subject to a variety of user-selectable parameters. 

2. [hivnetworkannotate](https://github.com/veg/hivclustering/wiki/hivnetworkannotate) : imports network node attributes (read from a `.csv`, `.tsv`, or `.json`) to a  `JSON` file output by `hivnetworkcsv`

3. [TNS](https://github.com/veg/hivclustering/wiki/TNS) : assuming sequences have associated sampling dates, computes the transmission network score (TNS), defined by [Little et al](https://www.ncbi.nlm.nih.gov/pubmed/24901437). Higher TNS is somewhat predictive of whether or not nodes in a network will accrue new edges in the future.