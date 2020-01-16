## Subcluster construction using `hivnetworkcsv`

> Sergei L Kosakovsky Pond (spond@temple.edu)

> 2020-01-09 v.1

Starting with version `0.6`, `hivnetworkcsv` allows the specification of multiple, comma-separated distance thresholds via the `-t` argument. This functionality permits inference of **nested** networks at different distance thresholds, i.e., the computation of subclusters of more "closely" linked nodes within clusters defined at a more permissive distance threshold. 

> ###This feature only affects runs where the JSON output (`-j` option) is specified.

To illustrate this feature, we will use the dataset `tests/subclusters/pirc.csv` which is a distance file constructed from `tests/subclusters/pirc.msa` -- sequences from [PMID 24901437](https://www.ncbi.nlm.nih.gov/pubmed/24901437).

Build the network, extracting date information from sequence names, and limiting the network to sequences from up to year 2009

```
hivnetworkcsv -j -i tests/subclusters/pirc.csv -t 0.02,0.01,0.005 -f regexp \
-p 0 "^([^|]+)\|([0-9]+)" -p 0 "^([^|]+)\|([0-9]+)" --before 20090101 \
-O tests/subclusters/pirc-2009.json
```

`stderr` output reports high level diagnostics.

```
At threshold 0.01 there were 77 subclusters
At threshold 0.005 there were 54 subclusters
```


This command builds the transmission network using the maximal thershold specified in the `-t` argument (0.02), and then creates further annotation of networks at 0.01 and 0.005 that are contained in the larger network. By construction, using the same input data, network at a **lower** distance threshold will be a subset of the netwotk inferred at a higher threshold. 

Additional annotation in the JSON file takes the following form.

1. A record in the `Settings` object<pre>
"Settings": {
    "threshold": 0.02,
    "edge-filtering": null,
    "contaminants": null,
    "compact_json": false,
    "created": "2020-01-09T22:10:33.833333+00:00",
    "additional thresholds": [
      0.01,
      0.005
    ]
  }
</pre>
2. A dictionary `Subclusters` which, for each threshold, reports clusters at that thershold. A cluster record is of the form <pre>
"ID" : [node1 index, node2 index, ... nodek index]
</pre> where `node_index` entries index into the `Node` array, also a part of the JSON output. For example, <pre>
"Subclusters": {
    "0.01": {
      "2.1": [
        24,
        25,
        26,
        27,
        28,
        29,
        30,
        31,
        32,
        33,
        34,
        35,
        36,
        37,
        38,
        39,
        40,
        41,
        42
      ],
      "19.1": [
        135,
        137
      ],
      "3.1": [
        43,
        44,
        45,
        46,
        47,
        48,
        49,
        50
      ],
      ...
 </pre>
 3. Individual node records in the `Node` array which belong to subclusters, are annotated as follows<pre>
 {
      "id": "KJ723207",
      "cluster": 9,
      "subcluster": [
        [
          0.01,
          "9.2"
        ],
        [
          0.005,
          "9.2"
        ]
      ]
      ...
    }
</pre> This node, in addition to belonging to cluster `9` at `t=0.02`, belongs to cluster `9.2` (a part of the larger cluster `9`) at both of the smaller thresholds.<pre>
{
      "id": "KJ723391",
      "cluster": 14,
      "baseline": [
        2002,
        1,
        1,
        0,
        0,
        0,
        1,
        1,
        -1
      ],
      "subcluster": [
        [
          0.01,
          "14.1"
        ]
      ]
    }
</pre> This node belongs to cluster `14` at `t=0.02`, cluster `14.1` at `t=0.01`, and no cluster at `t=0.005`

## Combining subcluster construction with incremental network building.

(See `incremental.md` for background on incremental construction).

When `-P` option is provided with `-t` and multiple thresholds, `hivnetworkcsv` will attempt to map subclusters from the current network to the subclusters from the previous network. The rule for matching subclusters is as follows.

If the current network has a subcluster, `C`, at distance threshold, `t`, **and** the previous network has a subcluster `Cp` at the same distance thereshold, then `C` will be matched with `Cp` *if and only if*

0. The subclusters each belong to the respective matched clusters.
1. The interesection of `C` and `Cp` is non-empty (i.e. **some** of the nodes from `C` are in `Cp`).
2. The intersection of `C` with any other previous network subcluster from the same cluster is empty, i.e. `C` does not overlap with any other existing subclusters. 
3. `Cp` does not overlap with any current network subclusters other that `C`

For example,

```
hivnetworkcsv -i tests/subclusters/pirc.csv -t 0.02,0.01,0.005 \
-f regexp -p 0 "^([^|]+)\|([0-9]+)" -p 0 "^([^|]+)\|([0-9]+)" \
-P tests/subclusters/pirc-2009.json -j -O tests/subclusters/pirc.json
```


