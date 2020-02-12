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
from math import log10, floor, sqrt, exp, log
from hivclustering import *
from functools import partial
import multiprocessing
import hppy as hy
from operator import itemgetter



run_settings = None
uds_settings = None

def settings():
    return run_settings

def ht_process_network_json (json):
    if 'trace_results' in json:
        json = json ['trace_results']
    if 'Settings' in json and 'compact_json' in json['Settings']:
        if json['Settings']['compact_json']:
            for key in ["Nodes","Edges"]:
                fields = list(json[key].keys())
                expanded = []
                for idx, f in enumerate (fields):
                    field_values = json[key][f]
                    if type (field_values) == dict and "values" in field_values:
                        field_values = [field_values["keys"][str(v)] for v in field_values["values"]]

                    for j,fv in enumerate(field_values):
                        if idx == 0:
                            expanded.append ({})
                        expanded[j][f] = fv

                json[key] = expanded

    return json

def ht_compress_network_json (network_info):
    def collect_keys (dict_set):
                unique_keys = set()
                for v in dict_set:
                    unique_keys.update (list (v.keys()))
                return unique_keys

    def compress_array (array):
        unique_values = {}
        try:
            for v in array:
                if v not in unique_values:
                    unique_values[v] = len (unique_values)

            #print (unique_values, file = sys.stderr)

            if len (unique_values) * 4 < len (array):
                lookup = {}
                for k, v in unique_values.items():
                    lookup[v] = k
                compact_array = {'keys' : lookup, 'values' : []}
                for a in array:
                    compact_array ['values'].append (unique_values[a])

                return compact_array

            return array

        except Exception as e:
            return array

    def convert_array_of_dicts (array, unique_keys):
        converted_set = {}
        null_by_key   = {}

        for k in unique_keys:
            converted_set[k] = []
            null_by_key[k] = 0


        for a in array:
            for k in unique_keys:
                if k in a:
                    if type (a[k]) is list:
                        converted_set[k].append (tuple (a[k]))
                    else:
                        converted_set[k].append (a[k])
                else:
                    converted_set[k].append (None)
                    null_by_key[k] += 1

        for k, nulls in null_by_key.items():
            '''
            print (k, nulls, len(array), file = sys.stderr)
            if nulls * 2 > len (array):
                sparse_set = {'sparse' : True, 'indices' : [], 'values' : []}
                for i, v in enumerate(converted_set[k]):
                    if v is not None:
                        sparse_set['indices'].append (i)
                        sparse_set['values'].append (v)
                sparse_set['values'] = compress_array(sparse_set['values'])
                converted_set[k] = sparse_set
            else:
            '''
            converted_set[k] = compress_array(converted_set[k])


        return converted_set

    network_info ["Edges"] = convert_array_of_dicts (network_info["Edges"], collect_keys (network_info["Edges"]))
    network_info ["Nodes"] = convert_array_of_dicts (network_info["Nodes"], collect_keys (network_info["Nodes"]))


def uds_attributes():
    return uds_settings

# Print iterations progress
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)

    https://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r', file = sys.stderr)
    # Print New Line on Complete
    if iteration == total:
        print(file = sys.stderr)

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
        print("Found %d clusters" % (len(clusters) - (1 if None in clusters else 0)), file=sys.stderr)
        cluster_sizes = sorted ([len(clusters[c]) for c in clusters if c is not None])
        print("Maximum cluster size = %d (second largest = %d) nodes" % (cluster_sizes[len(cluster_sizes)-1],cluster_sizes[len(cluster_sizes)-2]), file=sys.stderr)

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

    if not settings().skip_degrees:
        print("Fitting the degree distribution to various densities", file=sys.stderr)
        distro_fit = network.fit_degree_distribution()
        ci = distro_fit['rho_ci'][distro_fit['Best']]
        rho = distro_fit['rho'][distro_fit['Best']]
        bic = distro_fit['BIC'][distro_fit['Best']]
        rho = rho if rho is not None else 0.
        ci = ci if ci is not None else [0., 0.]
        if json_output:
            return_json['Degrees'] = {'Distribution': distro_fit['degrees'],
                                      'Model': distro_fit['Best'],
                                      'rho': rho,
                                      'rho CI': ci,
                                      'BIC': bic,
                                      'fitted': distro_fit['fitted'][distro_fit['Best']]}
        else:
            if (distro_fit['Best'] != "Negative Binomial"):
                ci = distro_fit['rho_ci'][distro_fit['Best']]
                rho = distro_fit['rho'][distro_fit['Best']]
                print("Best distribution is '%s' with rho = %g %s" %
                      (distro_fit['Best'], rho, ("[%g - %g]" % (ci[0], ci[1]))), file=sys.stderr)
            else:
                print("Best distribution is '%s'" % (distro_fit['Best']), file=sys.stderr)
    else:
        distro_fit = {'degrees' : network.get_degree_distribution()}

        if json_output:
            return_json['Degrees'] = distro_fit['degrees']

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
    for f in fn:
        fh = open(f)
        for line in fh:
            if line[0] == '>':
                yield line[1:].strip()

#-------------------------------------------------------------------------------

def compute_threshold_scores (full_records):

    def cluster_scaler (c):
        x = (cluster_max-c) / (cluster_max - cluster_min)
        return 1 - exp (- exp (-25*x + 3))

    def zscores (vector):
        mean = sum (vector) / len (vector)
        sigma = sqrt (sum ([(v-mean)**2 for v in vector]) / (len (vector)-1))
        zs = [(v-mean)/sigma for v in vector]
        zm = max (zs)
        return [z/zm for z in zs]

    records = []

    for line in full_records:
        records.append ([len (records), line[0], line[3], line[4] / max (1, line[5])])

    diffs  = []
    length = min (max (3, len(records)//100), 30)
    if length < 3:
        raise Exception ("Too few distance threshold datapoints (%d) to perform automatic threshold tuning" % len (full_records))

    cluster_min = min (records, key = lambda x: x[2])[2]
    cluster_max = max (records, key = lambda x: x[2])[2]

    for i, v in enumerate (records):
        if i >= length and i < len (records) - 1:
            trailing = sum ([k[3] for k in records [i - length : i + 1]])
            leading  = sum ([k[3] for k in records [i + 1: i + length + 1]])
            diffs.append ([i, trailing, leading, leading / trailing])

    diffs.sort (key = lambda r : r[3])
    zs = zscores ([d[3] for d in diffs])

    for i,d in enumerate(diffs):
        cs = cluster_scaler (records[d[0]][2])
        full_records[d[0]][6] = zs[i] + cs



#-------------------------------------------------------------------------------
def build_a_network(extra_arguments = None):

    random.seed()
    arguments = argparse.ArgumentParser(description='Construct a molecular transmission network.')

    arguments.add_argument('-i', '--input',   help='Input CSV file with inferred genetic links (or stdin if omitted). Can be specified multiple times for multiple input files (e.g. to include a reference database). Must be a CSV file with three columns: ID1,ID2,distance.', action = 'append')
    arguments.add_argument('-u', '--uds',   help='Input CSV file with UDS data. Must be a CSV file with three columns: ID1,ID2,distance.')
    arguments.add_argument('-d', '--dot',   help='Output DOT file for GraphViz (or stdout if omitted)')
    arguments.add_argument('-c', '--cluster', help='Output a CSV file with cluster assignments for each sequence')
    arguments.add_argument('-t', '--threshold', help='Only count edges where the distance is less than this threshold. Provide comma-separated values to compute subclusters if the output mode is JSON. If -t auto is specified, a heuristic is used to determine an ad hoc optimal threshold.')
    arguments.add_argument('-e', '--edi',   help='A .json file with clinical information')
    arguments.add_argument('-z', '--old_edi',   help='A .csv file with legacy EDI dates')
    arguments.add_argument('-f', '--format',   help='Sequence ID format. One of AEH (ID | sample_date | otherfiels default), LANL (e.g. B_HXB2_K03455_1983 : subtype_country_id_year -- could have more fields), regexp (match a regular expression, use the first group as the ID), or plain (treat as sequence ID only, no meta); one per input argument if specified', action = 'append')
    arguments.add_argument('-x', '--exclude',   help='Exclude any sequence which belongs to a cluster containing a "reference" strain, defined by the year of isolation. The value of this argument is an integer year (e.g. 1984) so that any sequence isolated in or before that year (e.g. <=1983) is considered to be a lab strain. This option makes sense for LANL or AEH data.')
    arguments.add_argument('-r', '--resistance',help='Load a JSON file with resistance annotation by sequence', type=argparse.FileType('r'))
    arguments.add_argument('-p', '--parser', help='The reg.exp pattern to split up sequence ids; only used if format is regexp; format is INDEX EXPRESSION (consumes two arguments)', required=False, type=str, action = 'append', nargs = 2)
    arguments.add_argument('-a', '--attributes',help='Load a CSV file with optional node attributes', type=argparse.FileType('r'))

    json_group = arguments.add_mutually_exclusive_group ();
    json_group.add_argument('-J', '--compact-json', dest = 'compact_json', help='Output the network report as a compact JSON object',required=False,  action='store_true', default=False)
    json_group.add_argument('-j', '--json', help='Output the network report as a JSON object',required=False,  action='store_true', default=False)

    arguments.add_argument('-o', '--singletons', help='Include singletons in JSON output',  action='store_true', default=False)
    arguments.add_argument('-k', '--filter', help='Only return clusters with ids listed by a newline separated supplied file. ', required=False)
    arguments.add_argument('-s', '--sequences', help='Provide the MSA with sequences which were used to make the distance file. Can be specified multiple times to include mutliple MSA files', required=False, action = 'append')
    arguments.add_argument('-n', '--edge-filtering', dest='edge_filtering', choices=['remove', 'report'], help='Compute edge support and mark edges for removal using sequence-based triangle tests (requires the -s argument) and either only report them or remove the edges before doing other analyses ', required=False)
    arguments.add_argument('-y', '--centralities', help='Output a CSV file with node centralities')
    arguments.add_argument('-l', '--edge-filter-cycles', dest = 'filter_cycles', help='Filter edges that are in cycles other than triangles', action='store_true')
    arguments.add_argument('--cycle-report-file', dest = 'cycle_report_filename', help='Prints cycle report to specified file', default = None, type = argparse.FileType('w'), required=False)
    arguments.add_argument('-g', '--triangles', help='Maximum number of triangles to consider in each filtering pass', type = int, default = 2**15)
    arguments.add_argument('-C', '--contaminants', help='Screen for contaminants by marking or removing sequences that cluster with any of the contaminant IDs (-F option) [default is not to screen]', choices=['report', 'remove'])
    arguments.add_argument('-F', '--contaminant-file', dest='contaminant_file',help='IDs of contaminant sequences', type=str)
    arguments.add_argument('-M', '--multiple-edges', dest='multiple_edges',help='Permit multiple edges (e.g. different dates) to link the same pair of nodes in the network [default is to choose the one with the shortest distance]', default=False, action='store_true')
    arguments.add_argument('-B', '--bridges',help='Report all bridges (edges whose removal would cause the graph to disconnect)', default=False, action='store_true')
    arguments.add_argument('--no-degree-fit', dest = "skip_degrees", help='Do not perform degree distribution fitting', default=False, action='store_true')
    arguments.add_argument('-X', '--extract',help='If provided, extract all the sequences ', type = int)
    arguments.add_argument('-O', '--output',help='Write the output file to', default = sys.stdout, type = argparse.FileType('w'))
    arguments.add_argument('-P', '--prior',help='When running in JSON output mode, provide a JSON file storing a previous (subset) version of the network for consistent cluster naming', required=False, type=argparse.FileType('r'))
    arguments.add_argument('-A', '--auto-profile', dest = 'auto_prof', help='If provided supercedes most other output and inference settings; will add edges from shortest to longest and report network statistics as a function of distance cutoff ', type = float)
    arguments.add_argument('--after', help='[assumes DATES are available] If provided (as YYYYMMDD) then only allow EDGES that connect nodes with dates at or AFTER this date', required=False, type = str)
    arguments.add_argument('--before', help='[assumes DATES are available] If provided (as YYYYMMDD) then only allow EDGES that connect nodes with dates at or BEFORE this date', required=False, type = str)
    arguments.add_argument('--import-attributes', dest = 'import_attr', help='Import node attributes from this JSON', required=False, type=argparse.FileType('r'))
    arguments.add_argument('--subcluster-annotation', dest = 'subcluster_annotation', help='As "dist" "field"". Use subcluster annotation for distance "dist" from node attribute "field"  ', required=False, nargs = 2)


    if extra_arguments:
        for a in extra_arguments:
            arguments.add_argument (*a["arg"], **a["kwarg"])

    global run_settings

    run_settings = arguments.parse_args()


    if run_settings.input == None:
        run_settings.input = [sys.stdin]
    else:
        try:
            run_settings.input = [open(file, 'r') for file in run_settings.input]
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

    formatter = []



    if run_settings.format is not None:
        regExpByIndex = {}

        if run_settings.parser is not None:
            for patterns in run_settings.parser :
                idx = int (patterns[0])
                if not idx in regExpByIndex:
                    regExpByIndex[idx] = []
                regExpByIndex[idx].append (patterns[1])

        for index, format_k in enumerate (run_settings.format):
            formats = {"AEH": parseAEH, "LANL": parseLANL, "plain": parsePlain, "regexp": parseRegExp(
                None if run_settings.parser is None  or index not in regExpByIndex else [re.compile (r) for r in regExpByIndex[index]])}
            try:
                formatter.append (formats[format_k])
            except KeyError:
                print("%s is not a valid setting for 'format' (must be in %s)" %
                      (run_settings.format, str(list(formats.keys()))), file=sys.stderr)
                raise

        if len (run_settings.format) != len (run_settings.input):
            raise Exception ("Must specify as many formatters as there are input files when at least one formatter is specified explicitly")
    else:
        formatter = [parseAEH for k in run_settings.input]

    if run_settings.exclude is not None:
        try:
            run_settings.exclude = datetime.datetime(int(run_settings.exclude), 12, 31)
        except ValueError:
            print("Invalid contaminant threshold year '%s'" % (run_settings.exclude), file=sys.stderr)
            raise

    run_settings.auto_threshold = False
    run_settings.additional_thresholds = None
    if run_settings.threshold is not None:
        if str (run_settings.threshold) == 'auto':
            run_settings.auto_threshold = True
            run_settings.threshold = None
        else:
            thresholds = [float (k) for k in run_settings.threshold.split (',')]
            run_settings.threshold = max (thresholds)
            if len (thresholds) > 1:
                run_settings.additional_thresholds = [k for k in thresholds if k != run_settings.threshold]
                run_settings.additional_thresholds.sort (reverse = True)

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

    if not run_settings.filter_cycles and run_settings.cycle_report_filename:
        raise ValueError('-l option is necessary to report cycles')

    network = transmission_network(multiple_edges=run_settings.multiple_edges)

    edge_filter_function = lambda edge : True

    if run_settings.before:
        run_settings.before = time.strptime(run_settings.before, '%Y%m%d')
        edge_filter_function = lambda edge : edge.check_exact_date (run_settings.before )

    if run_settings.after:
        run_settings.after = time.strptime(run_settings.after, '%Y%m%d')
        edge_filter_function = lambda edge, ef = edge_filter_function: ef (edge) and  edge.check_exact_date (run_settings.after, newer = True)

    if run_settings.auto_prof is not None or run_settings.auto_threshold:

        profile = []


        def network_report (threshold, network, max_clusters = [0]):
            clusters = network.retrieve_clusters(singletons=False)
            edges = len (network.edges)
            cl = sorted ([len (c) for c in clusters.values()], reverse = True)
            nnodes = sum (cl)
            profile.append ([threshold, sum (cl), edges, len (cl), cl[0] if len (cl) > 0 else 0, cl[1] if len (cl) > 1 else 0,0.])
            max_clusters[0] = max (max_clusters[0], len (cl))
            print('\rEvaluating distance threshold %8.5f %d %d' % (threshold, max_clusters[0], len (cl)), end = '\r', file = sys.stderr)

            #print ("%g\t%d\t%d\t%d\t%d\t%d\t%g" % (profile))
            sys.setrecursionlimit(max(sys.getrecursionlimit(), nnodes))
            return run_settings.auto_prof is not None or len (cl) == 0 or len (cl) > max_clusters[0] // 4

        network.read_from_csv_file_ordered(run_settings.input, network_report, formatter, run_settings.threshold if run_settings.threshold is not None else 1., 'BULK', run_settings.auto_prof if run_settings.auto_prof else 1e-5, filter = edge_filter_function)
        #print (profile, file = sys.stderr)
        compute_threshold_scores(profile)


        if run_settings.auto_prof is not None:
            print ("\t".join (["Threshold","Nodes","Edges","Clusters","LargestCluster","SecondLargestCluster","Score"]))
            for r in profile:
                print ("%g\t%d\t%d\t%d\t%d\t%d\t%g" % tuple(r))
            sys.exit (0)
        else:
            rec = [[k[0],k[-1]] for k in sorted (profile, key = lambda r : r[-1], reverse = True) if k [-1] >= 1.9]

            run_settings.threshold = None

            if len (rec) == 1:
                run_settings.threshold = rec[0][0]
            else:
                if len (rec) > 1:
                    suggested_span = max (rec, key = lambda x: x[0])[0] - min (rec, key = lambda x: x[0])[0]
                    mean_diff = sum ([k[1] - profile[i-1][1] for i,k in enumerate(profile[1:])]) / (len (profile)-1)
                    if (suggested_span / mean_diff < log (len (profile))):
                        run_settings.threshold = rec[0][0]

            if run_settings.threshold is None:
                if len (rec) == 0:
                    best_guess = sorted (profile, key = lambda r : r[-1], reverse = True)
                    print ('ERROR : Could not automatically determine a distance threshold; no sufficiently strong outlier, best guess %g (score %g)' % (best_guess[0][0], best_guess[0][-1]) , file = sys.stderr)
                else:
                    print ('ERROR : Multiple candidate thresholds: %s', ', '.join ([str (r[0]) for r in rec]), file = sys.stderr)
                sys.exit (1)
            else:
                print ("Selected distance threshold of % g" % run_settings.threshold, file = sys.stderr)
                to_delete = set ()
                for edge, v in network.edges.items():
                    if network.distances[edge] > run_settings.threshold:
                        to_delete.add (edge)
                for edge in to_delete:
                    del network.edges[edge]
                    del network.distances[edge]

    else:
        network.read_from_csv_file(run_settings.input, formatter, run_settings.threshold, 'BULK', filter = edge_filter_function)

    uds_settings = None

    if run_settings.uds:
        uds_settings = network.read_from_csv_file(run_settings.uds, formatter, run_settings.threshold, 'UDS', filter = edge_filter_function)

    sys.setrecursionlimit(max(sys.getrecursionlimit(), len (network.nodes)))

    if edi is not None:
        if old_edi:
            network.add_edi(edi)
        else:
            network.add_edi_json(edi)
        print("Added edi information to %d (of %d) nodes" %
              (len([k for k in network.nodes if k.edi is not None]), len (network.nodes)), file=sys.stderr)

        tabulate_edi = {}
        since_edi = [n.get_time_of_infection () for n in network.nodes]
        since_edi = ["Missing" if n is None else "Acute (<=90 days)" if n <= 90 else "Early (91-180 days)" if n <= 180 else "Chronic (>180 days)" for n in since_edi]
        for k in since_edi:
            if k not in tabulate_edi:
                tabulate_edi [k] = 1
            else:
                tabulate_edi [k] += 1

        for k in sorted (tabulate_edi.keys()):
            print("%s : %d" % (k, tabulate_edi[k]), file=sys.stderr)


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
            print ("Warning: Sequence ids referenced in input do not appear in source FASTA ids. Edge filtering will not be applied to all triangles.\n Missing ids in FASTA file: %s " %  ', '.join(missing_ids), file = sys.stderr)

        network.apply_attribute_filter('problematic', filter_out=True, do_clear=False)
        #if run_settings.filter:
        #    network.apply_id_filter(list=run_settings.filter, do_clear=False)

        network.compute_clusters () # this allocates nodes to clusters

        def generate_edges_by_cluster () :
            edges_by_clusters = {}
            current_edge_set = network.reduce_edge_set()
            for e in current_edge_set:
                cluster_id = e.p1.cluster_id
                if cluster_id not in edges_by_clusters:
                    edges_by_clusters[cluster_id] = [e]
                else:
                    edges_by_clusters[cluster_id].append (e)
            edges_by_clusters = [set(v) for c,v in edges_by_clusters.items() if len (v) >= 3]
            edges_by_clusters.sort (key = lambda x : len (x)) # smallest first
            return edges_by_clusters

        edges_by_clusters = generate_edges_by_cluster()


        # load sequence data

        hy_instance = hy.HyphyInterface()
        script_path = os.path.realpath(__file__)
        hbl_path = os.path.join(os.path.dirname(script_path), "data", "HBL", "ExtractSequences.bf")
        hy_instance.queuevar('_py_sequence_file', [os.path.abspath(s) for s in run_settings.sequences])
        hy_instance.runqueue(batchfile=hbl_path)

        all_referenced_sequences = set ()

        for cluster in edges_by_clusters:
            for e in cluster:
                all_referenced_sequences.update (e.sequences)

        referenced_sequence_data = {}

        for seq_id in all_referenced_sequences:
            referenced_sequence_data[seq_id] = hy_instance.getvar(seq_id, hy.HyphyInterface.STRING)

        if run_settings.extract is not None:
            sequence_set = set ()
            for anEdge in edges_by_clusters[run_settings.extract]:
                sequence_set.update (anEdge.sequences)
            for s in sequence_set:
                print (">%s\n%s\n" % (s, referenced_sequence_data[s]), file = run_settings.output)



        # partition edges into clusters

        def handle_a_cluster (edge_set, cluster_count, total_count):

            sys.stderr.write ('\r')
            sys.stderr.write ("Filtering a set of %d edges (%d/%d clusters) %s" % (len (edge_set), cluster_count, total_count, ' '*80))
            sys.stderr.flush ()

            edges_removed = 0
            #my_edge_set = edge_set
            maximum_number = run_settings.triangles

            supported_triangles = set ()

            for filtering_pass in range (8):
                edge_stats = network.test_edge_support(referenced_sequence_data, *network.find_all_simple_cycles(edge_set, maximum_number = maximum_number, ignore_this_set = supported_triangles), supported_cycles = supported_triangles)
                if not edge_stats:
                    break
                else:
                    edges_removed += edge_stats['removed edges']
                    #print("\tEdge filtering pass % d examined %d triangles, found %d poorly supported edges, and marked %d edges for removal" % (
                    #    filtering_pass, edge_stats['triangles'], edge_stats['unsupported edges'], edge_stats['removed edges']), file=sys.stderr)

                    sys.stderr.write ('\r')
                    sys.stderr.write ("Filtering a set of %d edges (%d/%d clusters). Pass %d, %d triangles, %d filtered edges" % (len (edge_set), cluster_count, total_count, filtering_pass, edge_stats['cycles'], edge_stats['removed edges']))
                    sys.stderr.flush ()
                    if edge_stats ['removed edges'] == 0:
                        break

                    maximum_number += run_settings.triangles
                    edge_set.difference_update (set ([edge for edge in edge_set if not edge.has_support()]))

            if run_settings.filter_cycles:
                maximum_number = run_settings.triangles
                supported_quads = set ()

                for filtering_pass in range (8):
                    simple_cycles = network.find_all_simple_cycles(edge_set, maximum_number = maximum_number, ignore_this_set = supported_quads, do_quads = True)

                    # If reporting cycle option set, pickle output to file
                    if run_settings.cycle_report_filename:
                        print(json.dumps({"cycles" : simple_cycles[0]}), file=run_settings.cycle_report_filename)

                    edge_stats = network.test_edge_support(referenced_sequence_data, *simple_cycles, supported_cycles = supported_quads, test_quads = True)

                    if not edge_stats:
                        break
                    else:
                        edges_removed += edge_stats['removed edges']
                        #print("\tEdge filtering pass % d examined %d triangles, found %d poorly supported edges, and marked %d edges for removal" % (
                        #    filtering_pass, edge_stats['triangles'], edge_stats['unsupported edges'], edge_stats['removed edges']), file=sys.stderr)

                        sys.stderr.write ('\r')
                        sys.stderr.write ("4-cycle filtering a set of %d edges (%d/%d clusters). Pass %d, %d 4-cycles, %d filtered edges" % (len (edge_set), cluster_count, total_count, filtering_pass, edge_stats['cycles'], edge_stats['removed edges']))
                        sys.stderr.flush ()

                        if edge_stats ['removed edges'] == 0:
                            break

                        maximum_number += run_settings.triangles
                        edge_set.difference_update (set ([edge for edge in edge_set if not edge.has_support()]))

            return edges_removed


        print ("Running edge filtering on %d clusters with 3 or more edges" % len (edges_by_clusters), file = sys.stderr)

        total_removed    = 0

        current_edge_set = set ()
        #require a minimim of X edges

        cluster_count    = 0
        max_edges_in_set = 0

        for edge_set in edges_by_clusters:
            current_edge_set.update (edge_set)
            max_edges_in_set = max (max_edges_in_set, len(edge_set))
            cluster_count += 1
            if max_edges_in_set >= 64:
                total_removed += handle_a_cluster (current_edge_set, cluster_count, len (edges_by_clusters))
                current_edge_set = set ()


        if len (current_edge_set) > 0:
            total_removed += handle_a_cluster (current_edge_set, cluster_count, len (edges_by_clusters))

        print ("\nEdge filtering identified %d edges for removal" % total_removed, file = sys.stderr)

        network.set_edge_visibility(edge_visibility) # restore edge visibility

        if run_settings.edge_filtering == 'remove':
            print("Edge filtering removed %d edges" % network.conditional_prune_edges(), file=sys.stderr)



    return network


