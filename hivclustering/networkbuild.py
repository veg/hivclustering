#!/usr/bin/env python3

import csv
import argparse
import operator
import sys
import datetime
import time
import random
import os.path
import json
import hppy as hy
import re
from math import log10, floor
from hivclustering import *
from functools import partial
import multiprocessing

run_settings = None
uds_settings = None


def settings():
    return run_settings


def uds_attributes():
    return uds_settings

#-------------------------------------------------------------------------------


def print_network_evolution(network, store_fitted=None, outdegree=False, distance=None, do_print=True, outfile=sys.stdout):
    byYear = []

    for year in range(2000, 2013):
        network.clear_filters()
        network.apply_date_filter(year, do_clear=True)
        if distance is not None:
            network.apply_distance_filter(distance, do_clear=False)
        network_stats = network.get_edge_node_count()
        network.compute_clusters()
        clusters = network.retrieve_clusters()
        if outdegree:
            distro_fit = network.fit_degree_distribution('outdegree')
        else:
            distro_fit = network.fit_degree_distribution()
        #print ("Best distribution is '%s' with rho = %g" % (distro_fit['Best'], 0.0 if distro_fit['rho'][distro_fit['Best']] is None else  distro_fit['rho'][distro_fit['Best']]), distro_fit['degrees'])
        if store_fitted is not None:
            store_fitted[year] = distro_fit['fitted']['Waring']
        byYear.append([year, network_stats['nodes'], network_stats['edges'], network_stats['total_sequences'], len(clusters), max(
            [len(clusters[c]) for c in clusters if c is not None]), distro_fit['rho']['Waring']] + distro_fit['rho_ci']['Waring'])

    #print (distro_fit)

    if do_print:
        print("\nYear,Nodes,Edges,Sequences,Clusters,MaxCluster,rho,rho_lower,rho_upper", file=outfile)
        for row in byYear:
            print(','.join([str(k) for k in row]), file=outfile)

#-------------------------------------------------------------------------------


def print_degree_distro(network, distro_fit, outfile=sys.stdout):
    print("\t".join(['degree', 'rawcount', 'rawpred', 'count', 'pred', 'ccount', 'cpred']), file=outfile)
    total = float(sum(distro_fit['degrees']))
    total1 = 0.
    total2 = 0.
    for k in range(0, len(distro_fit['degrees'])):
        vec = [str(p) for p in [k + 1, distro_fit['degrees'][k], distro_fit['fitted']['Waring'][k]
                                * total, distro_fit['degrees'][k] / total, distro_fit['fitted']['Waring'][k]]]
        vec.extend([0., 0.])
        total1 += distro_fit['degrees'][k] / total
        total2 += distro_fit['fitted']['Waring'][k]
        vec[5] = str(total1)
        vec[6] = str(total2)
        print("\t".join(vec))

    for dname, rho in distro_fit['rho'].items():
        print("%s : rho = %s, BIC = %s, p = %s" % (dname, 'N/A' if rho is None else "%5.2f" % (rho), 'N/A' if distro_fit["BIC"][
              dname] is None else "%7.2f" % distro_fit["BIC"][dname], 'N/A' if distro_fit["p"][dname] is None else "%4.2f" % (distro_fit["p"][dname])), file=outfile)


#-------------------------------------------------------------------------------
def describe_network(network, json_output=False, keep_singletons=False):
    network_stats = network.get_edge_node_count()
    if json_output:
        return_json = {'Network Summary': {'Edges': network_stats['edges'], 'Nodes': network_stats['nodes'],
                                           'Sequences used to make links': network_stats['total_sequences']},
                       'Multiple sequences': {'Subjects with': len(network_stats['multiple_dates']),
                                              'Followup, days': None if len(network_stats['multiple_dates']) == 0 else describe_vector([k[1] for k in network_stats['multiple_dates']])}
                       }

    else:
        print("%d edges on %d nodes" % (network_stats['edges'], network_stats['nodes']), file=sys.stderr)

    network.compute_clusters(keep_singletons)
    clusters = network.retrieve_clusters()
    #print (describe_vector([len(clusters[c]) for c in clusters]))

    if json_output:
        return_json['Network Summary']['Clusters'] = len(clusters)
        return_json['Cluster sizes'] = [len(clusters[c]) for c in clusters if c is not None]
    else:
        print("Found %d clusters" % len(clusters), file=sys.stderr)
        print("Maximum cluster size = %d nodes" % max([len(clusters[c])
                                                       for c in clusters if c is not None]), file=sys.stderr)

    if json_output:
        return_json['HIV Stages'] = {}

    if json_output:
        return_json['HIV Stages'] = network_stats['stages']
    else:
        for k in sorted (network_stats['stages'].keys()):
            print("%s : %d" % (k, network_stats['stages'][k]), file=sys.stderr)

    directed = 0
    reasons = {}
    for an_edge in network.reduce_edge_set():
        if an_edge.visible:
            if an_edge.compute_direction() is not None:
                directed += 1
            else:
                reason = an_edge.why_no_direction()
                if reason in reasons:
                    reasons[reason] += 1
                else:
                    reasons[reason] = 1

    if json_output:
        return_json['Directed Edges'] = {'Count': directed, 'Reasons for unresolved directions': reasons}
    else:
        print("%d directed edges" % directed, file=sys.stderr)
        print(reasons, file=sys.stderr)

    print("Fitting the degree distribution to various densities", file=sys.stderr)
    distro_fit = network.fit_degree_distribution()
    ci = distro_fit['rho_ci'][distro_fit['Best']]
    rho = distro_fit['rho'][distro_fit['Best']]
    rho = rho if rho is not None else 0.
    ci = ci if ci is not None else [0., 0.]
    if json_output:
        return_json['Degrees'] = {'Distribution': distro_fit['degrees'],
                                  'Model': distro_fit['Best'],
                                  'rho': rho,
                                  'rho CI': ci,
                                  'fitted': distro_fit['fitted'][distro_fit['Best']]}
    else:
        if (distro_fit['Best'] != "Negative Binomial"):
            ci = distro_fit['rho_ci'][distro_fit['Best']]
            rho = distro_fit['rho'][distro_fit['Best']]
            print("Best distribution is '%s' with rho = %g %s" %
                  (distro_fit['Best'], rho, ("[%g - %g]" % (ci[0], ci[1]))), file=sys.stderr)
        else:
            print("Best distribution is '%s'" % (distro_fit['Best']), file=sys.stderr)

    # find diffs in directed edges
    '''for anEdge in network.edges:
        if anEdge.visible:
            dir, diffr = anEdge.compute_direction (True)
            if dir is not None:
                print (diffr)
    '''
    if json_output:
        return return_json

    return distro_fit

#-------------------------------------------------------------------------------


def import_attributes(file, network):
    attribute_reader = csv.reader(file)
    header = next(attribute_reader)

    attribute_by_id = {}

    for line in attribute_reader:
        attribute_by_id[line[0]] = line[1]

    read_attributes = 0
    assigned = set()

    for a_node in network.nodes:
        if a_node.id in attribute_by_id:
            a_node.add_attribute(attribute_by_id[a_node.id])
            assigned.add(a_node.id)
            read_attributes += 1

    if read_attributes > 0:
        print('Loaded attribute information for %d/%d nodes' % (read_attributes, len(attribute_by_id)))
        print('Unassigned: ', set(attribute_by_id.keys()).difference(assigned))


#-------------------------------------------------------------------------------

def import_edi(file):
    edi_by_id = {}
    ediReader = csv.reader(file)
    header = next(ediReader)
    if len(header) != 14:
        raise Exception('Expected a .csv file with 14 columns as input')

    for line in ediReader:
        if len(line[1]):  # has PID
            id = line[1].replace('-', '')
        else:
            id = line[0]

        geno_date = None
        if len(line[2]):  # geno
            geno_date = time.strptime(line[2], '%m/%d/%Y')

        drug_date = None
        if len(line[4]):  # drugz
            drug_date = time.strptime(line[4], '%m/%d/%Y')

        edi_date = None
        stage = 'Chronic'

        if len(line[5]):  # disease stage
            stage = line[5]

        if len(line[6]):  # edi
            edi_date = time.strptime(line[6], '%m/%d/%Y')

        naive = False
        if line[3] == 'ARV Naive':
            naive = True

        if geno_date and edi_date:
            if edi_date > geno_date:
                # print time.mktime(edi_date) - time.mktime(geno_date)

                part1 = time.strftime("%m/%d", edi_date)
                part2 = time.strftime("%Y", geno_date)
                new_edi_date = time.strptime("/".join((part1, part2)), '%m/%d/%Y')
                #edi_date.tm_year = geno_date.tm_year
                if new_edi_date > geno_date:
                    continue
                else:
                    edi_date = new_edi_date

        viral_load = None
        if len(line[8]):  # vl
            viral_load = int(line[8])

        edi_by_id[id] = [geno_date, drug_date, stage, edi_date, viral_load, naive]
        #print (edi_by_id[id])
        # if (edi_date and drug_date and edi_date > drug_date):
        #	print "Fail %s" % id, edi_date, drug_date

    return edi_by_id

#-------------------------------------------------------------------------------


def import_edi_json(file):
    edi_by_id = json.load(file)
    for pid in edi_by_id:
        for key, value in edi_by_id[pid].items():
            if key == 'EDI':
                edi_by_id[pid]['EDI'] = time.strptime(edi_by_id[pid]['EDI'], '%Y-%m-%d')
            elif key == 'VL':
                for k in range(len(edi_by_id[pid]['VL'])):
                    edi_by_id[pid]['VL'][k][0] = tm_to_datetime(time.strptime(edi_by_id[pid]['VL'][k][0], '%Y-%m-%d'))
            elif key == 'ARV':
                edi_by_id[pid]['ARV'] = time.strptime(edi_by_id[pid]['ARV'], '%Y-%m-%d')
            else:
                edi_by_id[pid][key] = value

    return edi_by_id

#-------------------------------------------------------------------------------


def get_sequence_ids(fn):
    '''Expects newline separated file of node ids'''
    filter_list = set()
    with open(fn, 'r') as filter_file:
        reader = csv.reader(filter_file)
        for row in reader:
            filter_list.add(row[0])
        if not len(filter_list):
            pass
            #raise Exception('Empty file list')
    return filter_list


#-------------------------------------------------------------------------------

def get_fasta_ids(fn):
    fh = open(fn)
    for line in fh:
        if line[0] == '>':
            yield line[1:].strip()



#-------------------------------------------------------------------------------
def build_a_network():

    random.seed()
    arguments = argparse.ArgumentParser(description='Read filenames.')

    arguments.add_argument('-i', '--input',   help='Input CSV file with inferred genetic links (or stdin if omitted). Must be a CSV file with three columns: ID1,ID2,distance.')
    arguments.add_argument('-u', '--uds',   help='Input CSV file with UDS data. Must be a CSV file with three columns: ID1,ID2,distance.')
    arguments.add_argument('-d', '--dot',   help='Output DOT file for GraphViz (or stdout if omitted)')
    arguments.add_argument('-c', '--cluster', help='Output a CSV file with cluster assignments for each sequence')
    arguments.add_argument('-t', '--threshold', help='Only count edges where the distance is less than this threshold')
    arguments.add_argument('-e', '--edi',   help='A .json file with clinical information')
    arguments.add_argument('-z', '--old_edi',   help='A .csv file with legacy EDI dates')
    arguments.add_argument('-f', '--format',   help='Sequence ID format. One of AEH (ID | sample_date | otherfiels default), LANL (e.g. B_HXB2_K03455_1983 : subtype_country_id_year -- could have more fields), regexp (match a regular expression, use the first group as the ID), or plain (treat as sequence ID only, no meta)')
    arguments.add_argument('-x', '--exclude',   help='Exclude any sequence which belongs to a cluster containing a "reference" strain, defined by the year of isolation. The value of this argument is an integer year (e.g. 1984) so that any sequence isolated in or before that year (e.g. <=1983) is considered to be a lab strain. This option makes sense for LANL or AEH data.')
    arguments.add_argument('-r', '--resistance',help='Load a JSON file with resistance annotation by sequence', type=argparse.FileType('r'))
    arguments.add_argument('-p', '--parser', help='The reg.exp pattern to split up sequence ids; only used if format is regexp', required=False, type=str)
    arguments.add_argument('-a', '--attributes',help='Load a CSV file with optional node attributes', type=argparse.FileType('r'))
    arguments.add_argument('-j', '--json', help='Output the network report as a JSON object',required=False,  action='store_true', default=False)
    arguments.add_argument('-o', '--singletons', help='Include singletons in JSON output',required=False,  action='store_true', default=False)
    arguments.add_argument('-k', '--filter', help='Only return clusters with ids listed by a newline separated supplied file. ', required=False)
    arguments.add_argument('-s', '--sequences', help='Provide the MSA with sequences which were used to make the distance file. ', required=False)
    arguments.add_argument('-n', '--edge-filtering', dest='edge_filtering', choices=['remove', 'report'], help='Compute edge support and mark edges for removal using sequence-based triangle tests (requires the -s argument) and either only report them or remove the edges before doing other analyses ', required=False)
    arguments.add_argument('-y', '--centralities', help='Output a CSV file with node centralities')
    arguments.add_argument('-g', '--triangles', help='Maximum number of triangles to consider in each filtering pass', type = int, default = 2**16)
    arguments.add_argument('-C', '--contaminants', help='Screen for contaminants by marking or removing sequences that cluster with any of the contaminant IDs (-F option) [default is not to screen]', choices=['report', 'remove'])
    arguments.add_argument('-F', '--contaminant-file', dest='contaminant_file',help='IDs of contaminant sequences', type=str)
    arguments.add_argument('-M', '--multiple-edges', dest='multiple_edges',help='Permit multiple edges (e.g. different dates) to link the same pair of nodes in the network [default is to choose the one with the shortest distance]', default=False, action='store_true')

    global run_settings

    run_settings = arguments.parse_args()

    if run_settings.input == None:
        run_settings.input = sys.stdin
    else:
        try:
            run_settings.input = open(run_settings.input, 'r')
        except IOError:
            print("Failed to open '%s' for reading" % (run_settings.input), file=sys.stderr)
            raise

    if run_settings.dot is not None:
        try:
            run_settings.dot = open(run_settings.dot, 'w')
        except IOError:
            print("Failed to open '%s' for writing" % (run_settings.dot), file=sys.stderr)
            raise

    if run_settings.centralities is not None:
        try:
            run_settings.centralities = open(run_settings.centralities, 'w')
        except IOError:
            print("Failed to open '%s' for writing" % (run_settings.centralities), file=sys.stderr)
            raise

    edi = None
    old_edi = False

    if run_settings.edi is not None:
        try:
            run_settings.edi = open(run_settings.edi, 'r')
            edi = import_edi_json(run_settings.edi)
        except IOError:
            print("Failed to open '%s' for reading" % (run_settings.edi), file=sys.stderr)
            raise

    if edi is None and run_settings.old_edi is not None:
        try:
            run_settings.old_edi = open(run_settings.old_edi, 'r')
            edi = import_edi(run_settings.old_edi)
            old_edi = True
        except IOError:
            print("Failed to open '%s' for reading" % (run_settings.old_edi), file=sys.stderr)
            raise

    if run_settings.cluster is not None:
        try:
            run_settings.cluster = open(run_settings.cluster, 'w')
        except IOError:
            print("Failed to open '%s' for writing" % (run_settings.cluster), file=sys.stderr)
            raise

    formatter = parseAEH

    if run_settings.format is not None:
        formats = {"AEH": parseAEH, "LANL": parseLANL, "plain": parsePlain, "regexp": parseRegExp(
            None if run_settings.parser is None else re.compile(run_settings.parser))}
        try:
            formatter = formats[run_settings.format]
        except KeyError:
            print("%s is not a valid setting for 'format' (must be in %s)" %
                  (run_settings.format, str(list(formats.keys()))), file=sys.stderr)
            raise

    if run_settings.exclude is not None:
        try:
            run_settings.exclude = datetime.datetime(int(run_settings.exclude), 12, 31)
        except ValueError:
            print("Invalid contaminant threshold year '%s'" % (run_settings.exclude), file=sys.stderr)
            raise

    if run_settings.threshold is not None:
        run_settings.threshold = float(run_settings.threshold)

    if run_settings.uds is not None:
        try:
            run_settings.uds = open(run_settings.uds, 'r')
        except IOError:
            print("Failed to open '%s' for reading" % (run_settings.uds), file=sys.stderr)
            raise

    if len([k for k in [run_settings.contaminants, run_settings.contaminant_file] if k is None]) == 1:
        raise ValueError('Two arguments (-F and -S) are needed for contaminant screeening options')

    if len([k for k in [run_settings.edge_filtering, run_settings.sequences] if k is None]) == 1:
        raise ValueError('Two arguments (-n and -s) are needed for edge filtering options')

    network = transmission_network(multiple_edges=run_settings.multiple_edges)
    network.read_from_csv_file(run_settings.input, formatter, run_settings.threshold, 'BULK')

    uds_settings = None

    if run_settings.uds:
        uds_settings = network.read_from_csv_file(run_settings.uds, formatter, run_settings.threshold, 'UDS')

    sys.setrecursionlimit(max(sys.getrecursionlimit(), len (network.nodes)))

    if edi is not None:
        if old_edi:
            network.add_edi(edi)
        else:
            network.add_edi_json(edi)
        print("Added edi information to %d (of %d) nodes" %
              (len([k for k in network.nodes if k.edi is not None]), len (network.nodes)), file=sys.stderr)
        print("Added stage information to %d (of %d) nodes" %
              (len([k for k in network.nodes if k.stage is not None]), len (network.nodes)), file=sys.stderr)

    if run_settings.attributes is not None:
        import_attributes(run_settings.attributes, network)

    if run_settings.contaminant_file:
        run_settings.contaminant_file = get_sequence_ids(run_settings.contaminant_file)
        network.apply_cluster_membership_filter(run_settings.contaminant_file,
                                                filter_out=True, set_attribute='problematic')

        print("Marked %d nodes as being in the contaminant clusters" %
              len([n for n in network.nodes if n.has_attribute('problematic')]), file=sys.stderr)

        if run_settings.contaminants == 'remove':
            print("Contaminant linkage filtering removed %d edges" % network.conditional_prune_edges(
                condition=lambda x: x.p1.has_attribute('problematic') or x.p2.has_attribute('problematic')), file=sys.stderr)

    if run_settings.filter:
        run_settings.filter = get_sequence_ids(run_settings.filter)
        print("Included %d edges after applying node list filtering" %
              network.apply_cluster_membership_filter(run_settings.filter), file=sys.stderr)

    edge_visibility = network.get_edge_visibility()

    if run_settings.sequences and run_settings.edge_filtering:

        # Check that all sequences defined in distance file occur in source fasta file
        distance_ids = network.sequence_set_for_edge_filtering()
        source_fasta_ids = [id for id in get_fasta_ids(run_settings.sequences)]

        #print (distance_ids)

        if any(x not in source_fasta_ids for x in distance_ids):
            missing_ids = [x for x in distance_ids if x not in source_fasta_ids]
            raise Exception("Incorrect source file. Sequence ids referenced in input do not appear in source fasta ids. \n Missing ids in fasta file: %s " %  ', '.join(missing_ids))

        network.apply_attribute_filter('problematic', filter_out=True, do_clear=False)
        if run_settings.filter:
            network.apply_id_filter(list=run_settings.filter, do_clear=False)

        current_edge_set = network.reduce_edge_set()

        maximum_number = run_settings.triangles

        for filtering_pass in range (64):
            edge_stats = network.test_edge_support(os.path.abspath(
                run_settings.sequences), *network.find_all_triangles(current_edge_set, maximum_number = maximum_number))
            if not edge_stats or edge_stats['removed edges'] == 0:
                break
            else:
                print("Edge filtering pass % d examined %d triangles, found %d poorly supported edges, and marked %d edges for removal" % (
                    filtering_pass, edge_stats['triangles'], edge_stats['unsupported edges'], edge_stats['removed edges']), file=sys.stderr)

                maximum_number += run_settings.triangles
                current_edge_set = current_edge_set.difference (set ([edge for edge in current_edge_set if not edge.has_support()]))

        network.set_edge_visibility(edge_visibility)

        if edge_stats:
            print("Edge filtering examined %d triangles, found %d poorly supported edges, and marked %d edges for removal" % (
                edge_stats['triangles'], edge_stats['unsupported edges'], edge_stats['removed edges']), file=sys.stderr)
        else:
            print("Edge filtering examined %d triangles, found %d poorly supported edges, and marked %d edges for removal" % (
                0, 0, 0), file=sys.stderr)

        if run_settings.edge_filtering == 'remove':
            #print (len ([e for e in network.edge_iterator() if not e.has_support()]))
            print("Edge filtering removed %d edges" % network.conditional_prune_edges(), file=sys.stderr)
            # network.find_all_bridges()
    return network


