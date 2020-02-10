## Incremental network construction using `hivnetworkcsv`

> Sergei L Kosakovsky Pond (spond@temple.edu)

> 2020-01-09 v.1

For cases when a transmission network is built periodically from surveillance data, with a subsequent network (network `T`, for today) being (for the most part, see below) an expansion of a previous network (netwotk `Y`, for yesterday), `hivnetworkcsv` provides a **previous network** mode, which accomplishes the following, **assuming the `-j` (JSON output) option is selected**

1. Clusters in network `T` are matched, whenever practical, with clusters in network `Y`, so that _consistent_ cluster naming is maintained. 
2. Clusters in network `T` are **annotated** to explain how they compare to clusters in network `Y`. Annotation options for cluster `c` in network `T` are (see examples below). The annotation information is stored in the output JSON dictionary under the key `Cluster description`
	1. `existing` : cluster `c` is identical to a cluster in network `Y`	2. `expanded` : cluster `c` expands (adds node to) a cluster in network `Y`; this also includes the case of a contraction if some nodes were deleted from network `Y`.	3. `merged`  : cluster `c` merges two or more clusters from network `Y`, and possibly adds/removes nodes from them; example 
	4. `new`: cluster `c` consists of nodes that did not appear in clusters of network `Y`; this also includes the case when a cluster in network `Y` split into two more more smaller clusters in network `T`. 
3. Nodes in network `T` are **annotated** to explain how they relate to nodes in network `Y` (that these attributes are **not** mututally exclusive)
	1. `new_node` : this node does not appear in network `Y`
	2. `new_cluster` : this node belongs to a cluster in network `T` which has no match in cluster `Y`
	3. `moved_clusters`: this node was present in network `Y`, is present in network `T`, but moved to a different cluster
	4. *no attribute*: this node was present in network `Y`, is present in network `T`, and remains in the matched cluster
4. Edges in network `T` are **annotated** to explain how they relate to Edges in network `Y`
	1. `added-to-prior`: an edge in network `T` was not present in network `Y`. 

### Example

The `Y` (original) network is created by calling

```
hivnetworkcsv -t 0.015 -f plain -j -i tests/incremental/original.csv 
-O tests/incremental/original.json 
```

It consists of **5** clusters, with node names reflecting which cluster they belong to (e.g. `A1` belongs to cluster `1`)

**Figure 1** The existing network.

![Network `Y`](figures/incremental-1.png)

The `N` network is created in the incremental mode by calling.

```
hivnetworkcsv -t 0.015 -f plain -j -i tests/incremental/modified.csv \
-P tests/incremental/original.json -O tests/incremental/modified.json
```

> This example adds **and** removes nodes to the existing network


The following diagnostic messages are written to `stderr`

```
[WARNING] Removed 2 nodes from the previous network, E2, A1
Cluster 1 matches previous cluster 1
Cluster 2 matches previous cluster 1
[WARNING] Cluster 1 from the existing network is mapped to multiple new clusters / singletons (nodes E1, F1)
Cluster 3 MERGES 3, 4, None
Cluster 4 matches previous cluster 2
[WARNING] Incostistent singletons, nodes previously clustered became unclustered ({4})
Cluster 5 matches previous cluster 5
Cluster 6 is a new cluster
Added 4 edges compared to the prior network
```

**Figure 2** The updated network.

![Network `T`](figures/incremental-2.png)

**Figure 3** The changes from `Y` to `T`.

![Network `T`](figures/incremental-3.png)

1. Cluster `1` split in two clusters because node `A1` was deleted. The larger remaining cluster (3 nodes) inherited the cluster name (`1`), and the second remaining cluster (2 nodes, E1 and F1) was given a new name (`7`) and the nodes in that cluster were labeled as `moved`. In the output JSON this is annotated with 2. Cluster `2` lost one node (E2) due to deletion
3. Clusters `3` added a node (`E3`) and merged with cluster `4` through a link between the new node in cluster `4` (E4),  and cluster `3`. Cluster `3` also lost a node because the link between `B4` and `D4` was removed from network `T` (e.g., because the sequence for `D3` was revised). The resulting merged cluster inherited the ID of cluster `3`; because the existing clusters were of the same size and had to date information, the cluster with the lexicographically lowest node name (`A3`) inherited the name. Note, that cluster ID `4` has been retired, and will not appear in subsequent networks. 
4. Cluster `5` is unchanged.
5. Cluster `6` is new, i.e., none of its nodes were clustered in network `Y`.

 
### Additional information in JSON

Relevant sections of the output JSON are shown below

```
"Cluster description": {
    "1": {
      "type": "extended",
      "size": 3,
      "old_size": 6
    },
    "2": {
      "type": "extended",
      "size": 4,
      "old_size": 5
    },
    "3": {
      "type": "merged",
      "size": 9,
      "old_clusters": [
        3,
        4
      ],
      "sizes": [
        4,
        4
      ],
      "old_size": 8,
      "new_nodes": 2
    },
    "5": {
      "type": "existing"
    },
    "6": {
      "type": "new",
      "size": 2,
      "new_nodes": 2
    },
    "7": {
      "type": "new",
      "size": 2,
      "moved": 2
    }
  },
...

"Nodes" : [
...
	{
      "id": "A6",
      "cluster": 6,
      "attributes": [
        "new_node",
        "new_cluster"
      ]
	},
...
  {
      "id": "F1",
      "cluster": 7,
      "attributes": [
        "moved_clusters",
        "new_cluster"
      ]
    },
...

 {
      "id": "E4",
      "cluster": 3,
      "attributes": [
        "new_node"
      ]
    },
...
],

"Edges" : [
	...
	{
      "source": 7,
      "target": 8,
      "directed": false,
      "length": 0.01,
      "support": 0,
      "removed": false,
      "sequences": [
        "D3",
        "E3"
      ],
      "attributes": [
        "BULK",
        "added-to-prior"
      ]
    }
    ...
] 
 
```

