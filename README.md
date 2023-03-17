HIVClustering
-------------

A Python 3 library that makes infers molecular transmission networks from sequence data. A part of [HIV-TRACE](https://academic.oup.com/mbe/article/35/7/1812/4833215), available at http://hivtrace.org

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

```
usage: hivnetworkcsv [-h] [-i INPUT] [-u UDS] [-d DOT] [-c CLUSTER]
                     [-t THRESHOLD] [-e EDI] [-z OLD_EDI] [-f FORMAT]
                     [-x EXCLUDE] [-r RESISTANCE] [-p PARSER PARSER]
                     [-a ATTRIBUTES] [-J | -j] [-o {include,report}]
                     [-k FILTER] [-s SEQUENCES] [-n {remove,report}]
                     [-y CENTRALITIES] [-l]
                     [--cycle-report-file CYCLE_REPORT_FILENAME]
                     [-g TRIANGLES] [-C {report,remove}] [-F CONTAMINANT_FILE]
                     [-M] [-B] [--no-degree-fit] [-X EXTRACT] [-O OUTPUT]
                     [-P PRIOR] [-A AUTO_PROF] [--after AFTER]
                     [--before BEFORE] [--import-attributes IMPORT_ATTR]
                     [--subcluster-annotation SUBCLUSTER_ANNOTATION SUBCLUSTER_ANNOTATION]
                     [-q]

Construct a molecular transmission network.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input CSV file with inferred genetic links (or stdin
                        if omitted). Can be specified multiple times for
                        multiple input files (e.g. to include a reference
                        database). Must be a CSV file with three columns:
                        ID1,ID2,distance.
  -u UDS, --uds UDS     Input CSV file with UDS data. Must be a CSV file with
                        three columns: ID1,ID2,distance.
  -d DOT, --dot DOT     Output DOT file for GraphViz (or stdout if omitted)
  -c CLUSTER, --cluster CLUSTER
                        Output a CSV file with cluster assignments for each
                        sequence
  -t THRESHOLD, --threshold THRESHOLD
                        Only count edges where the distance is less than this
                        threshold. Provide comma-separated values to compute
                        subclusters if the output mode is JSON. If -t auto is
                        specified, a heuristic is used to determine an ad hoc
                        optimal threshold.
  -e EDI, --edi EDI     A .json file with clinical information
  -z OLD_EDI, --old_edi OLD_EDI
                        A .csv file with legacy EDI dates
  -f FORMAT, --format FORMAT
                        Sequence ID format. One of AEH (ID | sample_date |
                        otherfiels default), LANL (e.g. B_HXB2_K03455_1983 :
                        subtype_country_id_year -- could have more fields),
                        regexp (match a regular expression, use the first
                        group as the ID), or plain (treat as sequence ID only,
                        no meta); one per input argument if specified
  -x EXCLUDE, --exclude EXCLUDE
                        Exclude any sequence which belongs to a cluster
                        containing a "reference" strain, defined by the year
                        of isolation. The value of this argument is an integer
                        year (e.g. 1984) so that any sequence isolated in or
                        before that year (e.g. <=1983) is considered to be a
                        lab strain. This option makes sense for LANL or AEH
                        data.
  -r RESISTANCE, --resistance RESISTANCE
                        Load a JSON file with resistance annotation by
                        sequence
  -p PARSER PARSER, --parser PARSER PARSER
                        The reg.exp pattern to split up sequence ids; only
                        used if format is regexp; format is INDEX EXPRESSION
                        (consumes two arguments)
  -a ATTRIBUTES, --attributes ATTRIBUTES
                        Load a CSV file with optional node attributes
  -J, --compact-json    Output the network report as a compact JSON object
  -j, --json            Output the network report as a JSON object
  -o {include,report}, --singletons {include,report}
                        Include singletons in JSON output
  -k FILTER, --filter FILTER
                        Only return clusters with ids listed by a newline
                        separated supplied file.
  -s SEQUENCES, --sequences SEQUENCES
                        Provide the MSA with sequences which were used to make
                        the distance file. Can be specified multiple times to
                        include mutliple MSA files
  -n {remove,report}, --edge-filtering {remove,report}
                        Compute edge support and mark edges for removal using
                        sequence-based triangle tests (requires the -s
                        argument) and either only report them or remove the
                        edges before doing other analyses
  -y CENTRALITIES, --centralities CENTRALITIES
                        Output a CSV file with node centralities
  -l, --edge-filter-cycles
                        Filter edges that are in cycles other than triangles
  --cycle-report-file CYCLE_REPORT_FILENAME
                        Prints cycle report to specified file
  -g TRIANGLES, --triangles TRIANGLES
                        Maximum number of triangles to consider in each
                        filtering pass
  -C {report,remove}, --contaminants {report,remove}
                        Screen for contaminants by marking or removing
                        sequences that cluster with any of the contaminant IDs
                        (-F option) [default is not to screen]
  -F CONTAMINANT_FILE, --contaminant-file CONTAMINANT_FILE
                        IDs of contaminant sequences
  -M, --multiple-edges  Permit multiple edges (e.g. different dates) to link
                        the same pair of nodes in the network [default is to
                        choose the one with the shortest distance]
  -B, --bridges         Report all bridges (edges whose removal would cause
                        the graph to disconnect)
  --no-degree-fit       Do not perform degree distribution fitting
  -X EXTRACT, --extract EXTRACT
                        If provided, extract all the sequences
  -O OUTPUT, --output OUTPUT
                        Write the output file to
  -P PRIOR, --prior PRIOR
                        When running in JSON output mode, provide a JSON file
                        storing a previous (subset) version of the network for
                        consistent cluster naming
  -A AUTO_PROF, --auto-profile AUTO_PROF
                        If provided supercedes most other output and inference
                        settings; will add edges from shortest to longest and
                        report network statistics as a function of distance
                        cutoff
  --after AFTER         [assumes DATES are available] If provided (as
                        YYYYMMDD) then only allow EDGES that connect nodes
                        with dates at or AFTER this date
  --before BEFORE       [assumes DATES are available] If provided (as
                        YYYYMMDD) then only allow EDGES that connect nodes
                        with dates at or BEFORE this date
  --import-attributes IMPORT_ATTR
                        Import node attributes from this JSON
  --subcluster-annotation SUBCLUSTER_ANNOTATION SUBCLUSTER_ANNOTATION
                        As "dist" "field"". Use subcluster annotation for
                        distance "dist" from node attribute "field"
  -q, --quiet           Enable quiet mode


```
