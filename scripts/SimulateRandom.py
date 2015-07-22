#!/usr/bin/env python3.2

import csv, argparse, sys, datetime, time, random, os, math
#from scipy import stats
from hivclustering import *

N          = 648
replicates = 100
expected_stats = {'edges' : 540, 'nodes': 339, 'clusters' : 90}


def describe_network (network):
    network_stats = network.get_edge_node_count ()
    print ("%d edges on %d nodes" % (network_stats['edges'], network_stats['nodes']))
    network.compute_clusters()
    clusters = network.retrieve_clusters()
    print ("Found %d clusters" % len(clusters), file = sys.stderr)
    print ("Maximum cluster size = %d nodes" % max ([len (clusters[c]) for c in clusters if c is not None]), file = sys.stderr)

def crude_density_estimate (vector, min, max, step, value):
    dim = math.ceil((max-min)/step)
    density = [0.5 for k in range (dim)]
    
    for k in vector:
        bin = int ((k - min) // step)
        density[bin] += 1

    norm = sum (density)
    density = [k / norm  for k in density]
    
    return density [int ((value - min) // step)]

print ("\t".join(["Size","Bias","MeanEdges","MeanNodes","MeanClusters","P"]))

bias = 0

random.seed()

for int_bias in range (50,51,5):
    bias = int_bias/100.
    for size in range (1000,10001,250):
        edges    = []
        clusters = []
        nodes    = []
        
        for rep in range (replicates):
            sfn = transmission_network ()
            sfn = sfn.create_a_pref_attachment_network (size).sample_from_network (N,node_sampling_bias=bias)
            network_stats = sfn.get_edge_node_count ()
            sfn.compute_clusters()
            sfn_clusters = sfn.retrieve_clusters()
            edges += [network_stats['edges']]
            nodes += [network_stats['nodes']]
            clusters += [len (sfn_clusters)]
            #print (edges)
        
        prob = math.log(crude_density_estimate(edges,0,2000,10,expected_stats['edges'])*crude_density_estimate(nodes,0,650,5,expected_stats['nodes'])*crude_density_estimate(clusters,0,200,5,expected_stats['clusters']));
        
        print ("\t".join([str(k) for k in [size,bias,sum(edges)/replicates,sum(nodes)/replicates,sum(clusters)/replicates, prob]]))