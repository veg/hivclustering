#!/usr/bin/env python3

import json
import csv
import itertools
import os.path
import datetime

from .. import *
from ..networkbuild import *


def ensure_key (dict_object, key, default_value = {}):
    if key not in dict_object:
        dict_object[key] = default_value
    return dict_object[key]

def make_hiv_network():

    network = build_a_network()

    if settings().auto_prof:
        return network
        
        
    if settings().bridges:
        print ("Finding all non-trivial (not ending in a terminal node) bridges in the network", file = sys.stderr)
        network.find_all_bridges (attr = "bridge") ## this will labels the EDGES with "bridge" attribute

        bin_by_node = {}

        for e in network.reduce_edge_set():
            if e.has_attribute ("bridge") and e.p1.degree > 1 and e.p2.degree > 1 :
                if not settings().json:
                    print ("%s [degree %d] -- %s [degree %d]" % (e.p1.id, e.p1.degree, e.p2.id, e.p2.degree), file = settings().output)
                for n in [e.p1, e.p2]:
                    n.add_attribute ("bridge")
                    if n not in bin_by_node:
                        bin_by_node[n] = 1
                    else:
                        bin_by_node[n] += 1
            else:
                e.remove_attribute ("bridge")

        if not settings().json:
            if len (bin_by_node):
                for n, d in bin_by_node.items():
                    print ("Node %s (degree %d) is involved in %d bridges" % (n.id, n.degree, d), file = settings().output)
            else:
                print ("No non-trivial bridges found", file = sys.stderr)


    if settings().json or settings().compact_json:

        network_info = describe_network(network, True, settings().singletons)

        if settings().contaminant_file:
            network_info['Settings'] = {'threshold': settings().threshold,
                                        'edge-filtering': settings().edge_filtering,
                                        'contaminants': settings().contaminants,
                                        'contaminant-ids': list(settings().contaminant_file)
                                        }
        else:
            network_info['Settings'] = {'threshold': settings().threshold,
                                        'edge-filtering': settings().edge_filtering,
                                        'contaminants': settings().contaminants
                                        }

        if bool(settings().singletons):
            network_info['Settings']['singletons'] = True


        
        network_info['Settings']['compact_json'] = True if settings().compact_json else False
        network_info['Settings']['created'] = datetime.datetime.now (tz=datetime.timezone.utc).isoformat()

        nodes = []
        
        if settings().additional_thresholds:
            network_info ['Settings']['additional thresholds'] = settings().additional_thresholds
            network_info ["Subclusters"] = {}
            
            for e in network.reduce_edge_set():
                if not e.visible:
                    print (e)
            network.clear_adjacency()
            
            for n in network.nodes:
                n.original_cluster_id = n.cluster_id
            
            for t in settings().additional_thresholds:
                network.apply_distance_filter (t)
                network.compute_clusters()
                network_info ["Subclusters"][t] = network.retrieve_clusters (singletons = False)
                    
            network.clear_adjacency()
            for n in network.nodes:
                n.cluster_id = n.original_cluster_id
                del n.original_cluster_id
                
                
                
        cluster_summary_info = None
        prior_network_subclusters = None
        
        if settings().prior is not None:
            try:
                existing_nodes      = {} # node ID -> cluster ID
                prior_network       = ht_process_network_json(json.load (settings().prior))
                
                
                if "trace_results" in prior_network.keys():
                    prior_network_nodes = prior_network["trace_results"]["Nodes"]
                else:
                    prior_network_nodes = prior_network["Nodes"]

                existing_edges      = {} # node ID -> set of other node IDs
                prior_network_edges = prior_network["Edges"]
                new_clusters        = network.retrieve_clusters ()
                if 'Subclusters' in prior_network:
                    prior_network_subclusters = prior_network['Subclusters']
                else:
                    if settings().subcluster_annotation:
                        prior_network_subclusters = {}
                        subs = {}
                        field = settings().subcluster_annotation[1]
                        for i, n in enumerate (prior_network_nodes):
                            if field in n and "cluster" in n:
                                cluster_id = "%s.%s" % (n["cluster"], n[field])
                                ensure_key (subs, cluster_id, [])
                                subs[cluster_id].append (i)
                        prior_network_subclusters [settings().subcluster_annotation[0]] = subs
                    else:
                        prior_network_subclusters = None
                        
                max_id = 0
                cluster_sizes = {}

                old_cluster_sizes = {}
                new_cluster_sizes = {}
                node_to_prior_cluster = {}
                
                old_node_ids = set()


                for n in prior_network_nodes:
                    old_node_ids.add (n["id"])
                    if n["cluster"] is not None:
                        numeric_id = int (n["cluster"])
                        
                        node_to_prior_cluster [n["id"]] = numeric_id


                        if not numeric_id in old_cluster_sizes:
                            old_cluster_sizes[numeric_id] = 1
                        else:
                            old_cluster_sizes[numeric_id] += 1

                        existing_nodes[n["id"]] = numeric_id
                        max_id = max (max_id, numeric_id)
                        if n["cluster"] in cluster_sizes:
                            cluster_sizes[numeric_id] += 1
                        else:
                            cluster_sizes[numeric_id] = 1
                            
                            
                for e in prior_network_edges:
                    try:
                        node_pair = [prior_network_nodes[e['source']]["id"],prior_network_nodes[e['target']]["id"]]
                        for i, n in enumerate (node_pair):
                            if n not in existing_edges:
                                existing_edges[n] = set ()
                            existing_edges[n].add (node_pair [1-i])
                    except Exception as error:
                        if not settings().quiet: # new cluster
                            print ("[WARNING] Edge %s is not properly formed in the prior network" % str (e), file = sys.stderr)

                cluster_map = {} # // current cluster ID -> set of previous cluster ID
                
                # check if any of the old nodes got deleted
                
                new_node_ids = set ()

                for n in network.nodes:
                    new_node_ids.add (n.id)
                    if not n.cluster_id in new_cluster_sizes:
                        new_cluster_sizes[n.cluster_id] = 1
                    else:
                        new_cluster_sizes[n.cluster_id] += 1

                    if not n.cluster_id in cluster_map:
                        cluster_map[n.cluster_id] = set ()

                    if n.id in existing_nodes:
                        cluster_map[n.cluster_id].add (existing_nodes[n.id])
                    else:
                        cluster_map[n.cluster_id].add (None) # tag for a node not in previous clusters

                removed_nodes = old_node_ids - new_node_ids
                
                if len(removed_nodes):
                    ensure_key (network_info, "Notes",[])
                    network_info["Notes"].append ("Removed %d nodes from the previous network" % len(removed_nodes))
                    print ("[WARNING] %s, %s" % (network_info["Notes"][-1], ", ".join (list(removed_nodes))), file = sys.stderr)
                    
                mapped_to_existing_clusters = set ()
                used_existing_clusters = set ()

                cluster_remap = {}
                cluster_summary_info = {}
                '''
                    for each cluster id, record its properties
                '''
                discrepant_old_clusters = {}
                # this will categorize old clusters that didn't map to existing clusters via injection
                # id -> one or more of
                #    'multiple' : now maps to several new clusters [id1, id2]
                #    'disconnected': has nodes that have become singletons, T/F
                 

                def dump_cluster_members (id):
                    return ", ".join ([n.id for n in network.nodes if n.cluster_id == id]) 
                    
                def report_warning (message):
                    print ("[WARNING] %s" % message, file = sys.stderr)
                    
                def record_multiple_mapping (old_id, new_id):
                    multi_dict = ensure_key(ensure_key (discrepant_old_clusters, old_id, {}), 'multiple', set ())
                    multi_dict.add (new_id)
                    multi_dict.update (list(old_cluster_to_new_cluster[old_id]))
                    
                old_cluster_to_new_cluster = {}
                 
                for c_id, mapped_id in cluster_map.items():
                    if c_id is None: # singletons
                         if len (mapped_id) > 1 or None not in mapped_id:
                            report_warning ("Inconsistent singletons, nodes previously clustered became unclustered (%s)" % str (mapped_id))
                            for old_id in mapped_id:
                                ensure_key (discrepant_old_clusters, old_id)['disconnected'] = True
                                    
                    else:
                        if len (mapped_id) == 1:
                            if None in mapped_id :
                                if not settings().quiet: # new cluster
                                    print ("Cluster %d is a new cluster" % c_id,  file = sys.stderr)
                            else:
                                cluster_remap[c_id] = list(mapped_id)[0]
                                if not settings().quiet:
                                    print ("Cluster %d [%d nodes] matches previous cluster %d [%d nodes]" % (c_id, new_cluster_sizes [c_id], cluster_remap[c_id], old_cluster_sizes[cluster_remap[c_id]]),  file = sys.stderr)
                                cluster_summary_info [cluster_remap[c_id]] = {'type' : 'existing'}
                                if cluster_remap[c_id] in used_existing_clusters:
                                    record_multiple_mapping (cluster_remap[c_id], c_id)
                                    report_warning ("Cluster %d from the existing network is mapped to multiple new clusters / singletons (nodes %s)" % (cluster_remap[c_id], dump_cluster_members (c_id)))
                                else: # check to see if the older cluster shrunk
                                    if new_cluster_sizes [c_id] < old_cluster_sizes[cluster_remap[c_id]]:
                                        cluster_summary_info [cluster_remap[c_id]] = {'type' : 'extended', 'size' : new_cluster_sizes [c_id], 'old_size' : old_cluster_sizes[cluster_remap[c_id]]}
                                         
                                ensure_key (old_cluster_to_new_cluster, cluster_remap[c_id], set()).add (c_id)
                                used_existing_clusters.add (cluster_remap[c_id])
                                mapped_to_existing_clusters.add (c_id)
                        else:
                            if len (mapped_id) == 2 and None in mapped_id:
                                mapped_id.remove (None)
                                update_id = list(mapped_id)[0]
                                cluster_remap[c_id] = update_id
                                if not settings().quiet:
                                    print ("Cluster %d [%d nodes] extends previous cluster %d [%d nodes]" % (c_id, new_cluster_sizes [c_id], update_id, old_cluster_sizes[update_id]),  file = sys.stderr)
                                cluster_summary_info [update_id] = {'type' : 'extended', 'size' : new_cluster_sizes [c_id], 'old_size' : old_cluster_sizes[update_id]}
                                if cluster_remap[c_id] in used_existing_clusters:
                                    report_warning("Cluster %d from the existing network is mapped to multiple new clusters (nodes %s)" % (update_id, dump_cluster_members (c_id)))
                                    record_multiple_mapping (update_id, c_id)
                                    
                                ensure_key (old_cluster_to_new_cluster, update_id, set()).add (c_id)
                                used_existing_clusters.add (cluster_remap[c_id])
                                mapped_to_existing_clusters.add (c_id)
 
                            else:
                                if not settings().quiet:
                                    print ("Cluster %d MERGES %s" % (c_id, ", ".join ([str(k) for k in list(mapped_id)])), file = sys.stderr)
                                if None in mapped_id:
                                    mapped_id.remove (None)
                                    
                                if not mapped_id.isdisjoint (used_existing_clusters):
                                    report_warning ("One of the existing clusters being merged has already been mapped to a different new cluster (nodes: %s)" % dump_cluster_members (c_id))
                                    for old_id in mapped_id & used_existing_clusters:
                                        record_multiple_mapping (old_id, c_id)

                                
                                used_existing_clusters.update (mapped_id)
                                mapped_to_existing_clusters.add (c_id)
                                for old_id in mapped_id:
                                    ensure_key (old_cluster_to_new_cluster, old_id, set()).add (c_id)

                                mapped_id = sorted ([(k, cluster_sizes[k]) for k in mapped_id], key = lambda x: x[1], reverse = True)
                                largest_mapped_id = sorted ([k for k in mapped_id if k[1] == mapped_id[0][1]],  key = lambda x: x[0])
                                cluster_remap[c_id] = largest_mapped_id[0][0]
                                cluster_summary_info [largest_mapped_id[0][0]] = {'type' : 'merged',
                                                                          'size' : new_cluster_sizes [c_id],
                                                                          'old_clusters' : [k[0] for k in mapped_id],
                                                                          'sizes' : [k[1] for k in mapped_id],
                                                                          'old_size' : sum ([k[1] for k in mapped_id])}




                for d_id, d_items in discrepant_old_clusters.items():
                    '''
                        Strategies for resolving 'non-canonical clusters'
                            
                        1. If an old cluster is 'disconnected' only [no multiple]
                            This should mean that it is _completely_ disconnected, 
                            i.e. no connected components remain
                            Old code will handle it if a component of the cluster persists 
                                
                        2. If an old cluster maps to multiple new clusters, which themselves 
                           could map to multiple OLD clusters, then the following resolution algorithm applies
                           
                                -- If there are SOME NEW clusters that map to a SINGLE old cluster,
                                   then we choose one NEW cluster to inherit the ID of the old 
                                   cluster (using the standard sorting rules). The remaining
                                   new clusters will receive new IDs
                                   
                                   The single NEW cluster that "shrunk", or "shrunk" and added 
                                   new nodes, will be annotated as "part of old cluster" 
                                                                      
                                -- If ALL of the new clusters map to multiple OLD clusters, 
                                   then ALL of the  new clusters will receive new IDs
                           
                                   Nodes in NEW clusters that used to be in old clusters, 
                                   will be annotated as "migrated" [TODO, needs VIZ update]

                                   Clusters will be annotated as being parts of existing 
                                   clusters [TODO, needs VIZ update]
                    '''

                    if 'disconnected' in d_items:
                        if not 'multiple' in d_items:
                            #print ("Disconnected and multiple")
                            continue
                    else:
                        if 'multiple' in d_items:  
                            # find to see if there are any NEW clusters that map ONLY to this old cluster (and maybe Nones)
                            new_mapping_only_to_old = [c_id for c_id in d_items['multiple'] if len (cluster_map[c_id]) == 1 or len (cluster_map[c_id]) == 2 and None in cluster_map[c_id]]
                            if len (new_mapping_only_to_old): # select a cluster to inherit the old ID
                                if len (new_mapping_only_to_old) == 1:
                                    inheritor = list (new_mapping_only_to_old)[0]
                                else:
                                    candidate_clusters = {}
                                    for ncid in new_mapping_only_to_old:
                                        candidate_clusters[ncid] = new_clusters[ncid]
                                    inheritor = network.sort_clusters(singletons = False, precomputed_clusters = candidate_clusters, set_cluster_id = lambda node, value: node)[1][0].cluster_id
                                    
                                
                                # update the record for the inheritor (which has shrunk)
                                cluster_summary_info [d_id] = {'type' : 'extended', 'size' : new_cluster_sizes [inheritor], 'old_size' : old_cluster_sizes[d_id]}
                                              
                                              
                                leftovers = set ()
                                leftovers.update (d_items['multiple']) 
                                # these are NEW clusters that are mapped to THIS old cluster AND some other clusters
                                # they need to be given a new ID               
                                leftovers.remove (inheritor)                                 
                                for nid in new_mapping_only_to_old:
                                    if nid != inheritor:  
                                        leftovers.remove (nid)
                                        if nid in cluster_remap:
                                            del cluster_remap[nid]
                                        for n in new_clusters[nid]:
                                            if n.id in old_node_ids:
                                                n.add_attribute ("moved_clusters")
                                                
                                for nid in leftovers:
                                    if nid in cluster_remap:
                                        del cluster_remap[nid]
                                    for n in new_clusters[nid]:
                                        if n.id in old_node_ids:
                                            n.add_attribute ("moved_clusters")
                                
                            else:
                                for nid in d_items['multiple']:
                                    if nid in cluster_remap:
                                        del cluster_remap[nid]
                                    for n in new_clusters[nid]:
                                        if n.id in old_node_ids:
                                            n.add_attribute ("moved_clusters")
                                        
                            

 
                    #print (d_id, d_items, file = sys.stderr)
                    #for c_id in d_items['multiple']:
                    #    print (cluster_map[c_id], file = sys.stderr)
                        
                        
                        
        
                
                for n in network.nodes:
                    if not n.id in existing_nodes:
                        n.add_attribute ("new_node")
                        node_to_prior_cluster[n.id] = None
                    else:
                        n.remove_attribute ("new_node")
                        
                    n.remove_attribute ("new_cluster")

                    if n.cluster_id in cluster_remap:
                        n.cluster_id  = cluster_remap [n.cluster_id]
                    else:
                        if n.cluster_id is not None: # new clusters go here
                            n.cluster_id  = -n.cluster_id
                            n.add_attribute ("new_cluster")

                #print (mapped_to_existing_clusters, "\n\n", cluster_remap, file = sys.stderr)
                cluster_info = network.sort_clusters(filter = lambda cluster_id ,cluster_data : cluster_id < 0, start_id = max_id + 1)

                for n in network.nodes:
                    if n.cluster_id is not None:
                        if n.has_attribute ("new_cluster"):
                            if n.cluster_id not in cluster_summary_info:
                                cluster_summary_info[n.cluster_id] = {'type' : 'new', 'size' : 1}
                            else:
                                #print (n.cluster_id, max_id, cluster_summary_info[n.cluster_id], file = sys.stderr)
                                cluster_summary_info[n.cluster_id]['size'] += 1
                            
                        if n.has_attribute ("new_node"):
                            ensure_key (cluster_summary_info[n.cluster_id], 'new_nodes', 0)
                            cluster_summary_info[n.cluster_id]['new_nodes'] += 1
                        
                        if n.has_attribute ("moved_clusters"):
                             if n.cluster_id not in cluster_summary_info:
                                raise Exception ("Internal error in cluster annotation")
                             else:
                                ensure_key (cluster_summary_info[n.cluster_id], 'moved', 0)
                                cluster_summary_info[n.cluster_id]['moved'] += 1
                            
                       

                added_edges = 0
    
                for e in network.edges:
                    p1p2 = False
                    p2p1 = False
                    if e.p1.id in existing_edges:
                        if not e.p2.id in existing_edges[e.p1.id]:
                            e.update_attributes ("added-to-prior")
                        else:  
                            p1p2 = True
                    if e.p2.id in existing_edges:
                        if not e.p1.id in existing_edges[e.p2.id]:
                            e.update_attributes ("added-to-prior")
                        else:  
                            p2p1 = True
                            
                    if p1p2 != p2p1: 
                        print ("Edge ", e, " is not consistent with prior network", file = sys.stderr)
                    else:
                        if not p1p2:
                            added_edges += 1
                            
                print ("Added %d edges compared to the prior network" % added_edges, file = sys.stderr)
                        
            except Exception as e:
                print ("Error with prior network processing: %s" % str (e), file = sys.stderr)
                raise
                sys.exit (1)

        else:
            cluster_info = network.sort_clusters()

        node_idx = {}

        for idx, cluster in cluster_info.items():
            #print (idx, cluster)
            for n in cluster:
                if idx is not None:
                    node_info = {'id' : n.id, 'cluster' : idx}
                    if n.attributes and len (n.attributes):
                        node_info['attributes'] = list(n.attributes)
                    if n.get_edi():
                        node_info['edi'] = n.get_edi()
                    bd = n.get_baseline_date(True)
                    if bd:
                        node_info['baseline'] = bd
                    n.cluster_id = idx # this is for handling subclusters later
                          
                    nodes.append(node_info)
                    node_idx[n] = len(nodes) - 1



        edges = []
        for e in network.reduce_edge_set():
            if e.visible:
                edge_source = e.compute_direction()
                if edge_source is not None:
                    src = node_idx[edge_source]
                    rcp = node_idx[e.p2 if edge_source != e.p2 else e.p1]
                    directed = True
                else:
                    src = node_idx[e.p1]
                    rcp = node_idx[e.p2]
                    directed = False
                edges.append({'source': src, 'target': rcp, 'directed': directed, 'length': network.distances[
                             e], 'support': e.edge_reject_p, 'removed': not e.has_support(), 'sequences': e.sequences, 'attributes': list(e.attribute)})

        network_info['Nodes'] = nodes
        # do NOT sort the nodes, otherwise referencing them by index from edges will be broken.
        network_info['Edges'] = sorted(edges, key=lambda edge: ''.join([sequence for sequence in edge['sequences']]))
        # OK to sort edges, because their order is not used downstream
        if cluster_summary_info is not None:
            network_info['Cluster description'] = cluster_summary_info
            
        
        if 'Subclusters' in network_info:
        
            
            if prior_network_subclusters:
               prior_network_subclusters_by_id = {}
               for distance, cluster in prior_network_subclusters.items():
                   for id, members in cluster.items():
                        cluster_id = id.split (".")[0]
                        subscluster_id = id.split (".")[1]
                        ensure_key (prior_network_subclusters_by_id, cluster_id, {})
                        ensure_key (prior_network_subclusters_by_id[cluster_id], distance, {})
                        prior_network_subclusters_by_id[cluster_id][distance][subscluster_id] = frozenset ([prior_network_nodes[k]["id"] for k in members])
            else:
                prior_network_subclusters_by_id = None       
                   
        
            reindexed_subclusters = {}
            
            for t, subcluster_set in network_info["Subclusters"].items():
                cluster_to_subcluster_idx = {}
                reindexed_subclusters_t   = {}
                for id, node_set in subcluster_set.items():
                    if subcluster_set[id][0].cluster_id not in cluster_to_subcluster_idx:
                        cluster_to_subcluster_idx[ subcluster_set[id][0].cluster_id] = {}
                    cluster_to_subcluster_idx[subcluster_set[id][0].cluster_id][id] = node_set
                     
                for cluster_id, subclusters in cluster_to_subcluster_idx.items(): 
                    start_id = 1
                    cstr = str (cluster_id)
                    tstr = str (t)
                    matched_subclusters = {}
                    reverse_matched = {}
                    cluster_filter = None
                    if prior_network_subclusters_by_id and cstr in prior_network_subclusters_by_id and tstr in prior_network_subclusters_by_id[cstr]:
                        start_id = 1 + len (prior_network_subclusters_by_id[cstr][tstr])                    
                        for sub_id,subcluster_nodes in subclusters.items():
                            my_nodes = frozenset ([n.id for n in subcluster_nodes])
                            for pid, previous_set in prior_network_subclusters_by_id[cstr][tstr].items():
                                if not my_nodes.isdisjoint (previous_set):
                                    ensure_key (matched_subclusters, sub_id, [])
                                    matched_subclusters[sub_id].append (pid)
                                    ensure_key (reverse_matched, pid, [])
                                    reverse_matched[pid].append (sub_id)
                                    
                            cluster_filter = lambda id,nodes: id not in matched_subclusters or len(matched_subclusters[id]) != 1 or len (reverse_matched[matched_subclusters[id][0]]) !=1
                                    
                        for sub_id, matches in matched_subclusters.items():
                            if len (matches) == 1 and len (reverse_matched[matches[0]]) == 1:
                                for n in subclusters[sub_id]:
                                    setattr(n, "subcluster_id", int (matches[0]))
                        
                        
                    for sub_id,subcluster_nodes in network.sort_clusters(singletons = False, precomputed_clusters = subclusters, filter = cluster_filter, start_id = start_id, set_cluster_id = lambda node, value: setattr(node, "subcluster_id", value)).items():
                        subcluster_id = subcluster_nodes[0].subcluster_id
                        combined_id = "%d.%d" % (cluster_id, subcluster_id)
                        reindexed_subclusters_t[combined_id] = [node_idx[node] for node in subcluster_nodes]
                        for node in subcluster_nodes: 
                            node_for_json = nodes[node_idx[node]]
                            if not 'subcluster' in node_for_json:
                                node_for_json['subcluster'] = []
                            node_for_json['subcluster'].append ((t,combined_id))
                            
                             
                        
                reindexed_subclusters[t]    = reindexed_subclusters_t
                
            network_info['Subclusters'] = reindexed_subclusters
            
            for t, subcluster_set in network_info["Subclusters"].items():
                print ("At threshold %g there were %d subclusters" % (t, len(subcluster_set)), file = sys.stderr)
                
        if settings().import_attr:
            import_data = json.load (settings().import_attr)
            if 'trace_results' in import_data:
                import_data = import_data['trace_results']
            
            if 'patient_attribute_schema' in import_data and 'Nodes' in import_data:
                network_info['patient_attribute_schema'] = import_data['patient_attribute_schema']
                id_to_idx = {}
                for i,n in enumerate(network_info["Nodes"]):
                    id_to_idx[n['id']] = i
                for n in import_data['Nodes']:
                    if n['id'] in id_to_idx and 'patient_attributes' in n:
                        network_info["Nodes"][id_to_idx[n['id']]]['patient_attributes'] = n['patient_attributes']
                
        if settings().compact_json:
            ht_compress_network_json (network_info)
            
  
        print(json.dumps(network_info, separators=(',', ':')), file = settings().output)

    else:
        describe_network(network)

    if settings().dot:
        network.generate_dot(settings().dot)

    if settings().cluster:
        network.write_clusters(settings().cluster, node_to_prior_cluster if settings().prior is not None and settings().json else None)

    if settings().centralities:
        network.write_centralities(settings().centralities)

    return network

if __name__ == '__main__':
    __spec__ = None
    make_hiv_network()

