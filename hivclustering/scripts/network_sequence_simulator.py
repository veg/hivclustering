
import json
import csv
import os.path
import hashlib
import subprocess
from math import sqrt

def positive_integer (value):
    ivalue = int(value)
    if ivalue <= 0:
         raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
    return ivalue

def float01 (value):
    fvalue = float(value)
    if fvalue < 0. or fvalue > 1.:
         raise argparse.ArgumentTypeError("%s is an invalid rate " % value)
    return fvalue

def nn_float (value):
    fvalue = float(value)
    if fvalue < 0.:
         raise argparse.ArgumentTypeError("%s is an invalid rate " % value)
    return fvalue


from Bio import SeqIO

from hivclustering import *
from hivclustering.networkbuild import *

def mcc (m):
    return (m[1][1]*m[0][0] - m[1][0]*m[0][1])/sqrt ((m[1][1] + m[0][1])*(m[1][1] + m[1][0])*(m[0][0] + m[0][1])*(m[0][0] + m[1][0]))

def run_tn93 (in_path, out_path, threshold = 0.015):

    try:
        subprocess.check_call (['/usr/local/bin/tn93', '-q', '-t', str(threshold), '-o', out_path, in_path]) 
    except subprocess.CalledProcessError as e:
        print ('ERROR: tn93 call failed',e,file = sys.stderr)
    
    return None

if __name__=='__main__':
    random.seed()
    arguments = argparse.ArgumentParser(description='Read filenames.')
    arguments.add_argument('-s', '--sequences', help = 'Provide the MSA with sequences which were used to make the distance file. ', required = True)
    arguments.add_argument('-f', '--fasta', help = 'Write simulated network to. ', required = True)
    arguments.add_argument('-t', '--tn93', help = 'Write the CSV file to. ', required = True, type = str)
    arguments.add_argument('-n', '--size', help = 'Number of sequences to simulate ', required = False, type = positive_integer, default = 500)
    arguments.add_argument('-d', '--days', help = 'Mean number of days before transmissions', required = False, type = nn_float, default = 7)
    arguments.add_argument('-r', '--rate', help = 'Expected subs/site/year', required = False, type = float01, default = 0.0007)
    arguments.add_argument('-l', '--lineages', help = 'Starting lineages', required = False, type = positive_integer, default = 10)
    arguments.add_argument('-p', '--split', help = 'Initiate lineages at this rate', required = False, type = float01, default = 0.02)
    arguments.add_argument('-m', '--random', help = 'Make random attachments at this rate', required = False, type = nn_float, default = 0.2)
    arguments.add_argument('-u', '--subset', help = 'Subsample this many sequences', required = False, type = positive_integer)
    arguments.add_argument('-b', '--bias', help = 'Bias subsamples to link to already sampled nodes', required = False, type = nn_float, default = 0.0)
    arguments.add_argument('-y', '--sampling', help = 'Mean number of days before sampling', required = False, type = positive_integer, default = 30)
    arguments.add_argument('-x', '--burst', help = 'Mean (poisson) number of individuals infected at a given date', required = False, type = nn_float, default = None)
    arguments.add_argument('-e', '--replicates', help = 'Simulate this many replicates', required = False, type = positive_integer, default = 1)
    arguments.add_argument('-T', '--threshold', help = 'Distance threshold for connecting edges.', required = False, type = float, default = 0.015)
    arguments.add_argument('-F', '--edge-filtering', dest = 'edge_filtering', help = 'Apply edge filtering (false by default).', action = 'store_true', default = False)

    
    settings = arguments.parse_args()
    
    ambig_remover = str.maketrans ("RYSWKMBDHVN-",
                                   "ACCAGACAAAAA")
    
    with open (settings.sequences, "r") as fh:
        sequences =  [str(record.seq).translate (ambig_remover)[0:1212] for record in SeqIO.parse (fh, "fasta")]
        
        

    print ("Read %d sequences " % len(sequences))
        
    kendall_p_values = []
    ppv              = []
        
    sys.setrecursionlimit(max(sys.getrecursionlimit(),settings.size))

    for replicate in range (settings.replicates): 
        print ("##### REPLICATE %d #####" % (replicate+1))  
        random_network = transmission_network ()

        start_nodes = random_network.create_a_pref_attachment_network (network_size = settings.size, start_with = settings.lineages, random_attachment = settings.random, start_new_tree = settings.split, start_date = datetime.datetime (1996,1,1), tick_rate = settings.days, poisson_mean = settings.burst)
        #def sample_from_network (self, how_many_nodes = 100, how_many_edges = None, node_sampling_bias = 0.0):

        if settings.subset is not None: 
            subset_network = random_network.sample_from_network (settings.subset,node_sampling_bias = settings.bias)
        else:
            subset_network = random_network
        
        json1 = describe_network (random_network, json_output = True)
        print ("Nodes %d, edges %d, clusters %d, distro %s, rho %g" % 
                (json1['Network Summary']['Nodes'],json1['Network Summary']['Edges'],json1['Network Summary']['Clusters'],json1['Degrees']['Model'], json1['Degrees']['rho']))
                
        given_degrees = random_network.get_node_degree_list ()
        
        for e in random_network.edges:
            print (e)
        
        #sys.exit (0)
    
        #describe_network (random_network)
        print ("Sampling sequences...", file = sys.stderr)
        seqs = random.sample (sequences, len (start_nodes))
    
        random_network.simulate_sequence_evolution (start_nodes, seqs, settings.rate, settings.sampling)
        with open (settings.fasta, 'w') as fh:
            random_network.dump_as_fasta (fh, add_dates = False, filter_on_set = subset_network.nodes if settings.subset is not None else None)
    
        #def sample_from_network (self, how_many_nodes = 100, how_many_edges = None, node_sampling_bias = 0.0):
     
    
        #print ("\n\n", given_degrees)
        run_tn93 (settings.fasta, settings.tn93, settings.threshold)
    
        recovered_network = transmission_network ()
        with open (settings.tn93) as fh:
            recovered_network.read_from_csv_file (fh, parsePlain, settings.threshold, 'BULK')
            
        if settings.edge_filtering:
        
            edge_stats = recovered_network.test_edge_support(os.path.abspath(
                settings.fasta), *recovered_network.find_all_triangles(recovered_network.reduce_edge_set()))

            if edge_stats:
                print("Edge filtering examined %d triangles, found %d poorly supported edges, and marked %d edges for removal" % (
                    edge_stats['triangles'], edge_stats['unsupported edges'], edge_stats['removed edges']), file=sys.stderr)
            else:
                print("Edge filtering examined %d triangles, found %d poorly supported edges, and marked %d edges for removal" % (
                    0, 0, 0), file=sys.stderr)

            with open ("/tmp/recovered.dot", "w") as fh:
                recovered_network.generate_dot (fh, attribute_color = lambda x: "red" if x.edge_reject_p > 0.05 else None)   
            
            print("Edge filtering removed %d edges" % recovered_network.conditional_prune_edges(), file=sys.stderr)
        
        json2 = describe_network (recovered_network, json_output = True)
        print ("Nodes %d, edges %d, clusters %d, distro %s, rho %g" % 
                (json2['Network Summary']['Nodes'],json2['Network Summary']['Edges'],json2['Network Summary']['Clusters'],json2['Degrees']['Model'], json2['Degrees']['rho']))
                

        for e in recovered_network.edges:
            print (e, " [X] " if e in random_network.edges else "")
                
        with open ("/tmp/random.dot", "w") as fh:
            random_network.generate_dot (fh)   

        with open ("/tmp/filtered.dot", "w") as fh:
            recovered_network.generate_dot (fh)   
            


