

import datetime, time, random, itertools, operator, re, sys
from math import log
from copy import copy, deepcopy
from bisect import bisect_left
from operator import itemgetter
import hppy as hy
import os
import csv
import multiprocessing
from functools import partial

__all__ = ['edge', 'patient', 'transmission_network', 'parseAEH', 'parseLANL', 'parsePlain', 'parseRegExp', 'describe_vector', 'tm_to_datetime', 'datetime_to_tm']
#-------------------------------------------------------------------------------


def parseAEH (str):
    try:
        bits = str.rstrip().split ('|')
        if len (bits) < 2:
            raise Exception ('Improperly formatted AEH header (need at least "ID|Sample date in mmddyyyy format": %s' % str)

        patient_description = {}
        patient_description ['id']   =  bits[0]
        patient_description ['date'] = time.strptime (bits[1],'%m%d%Y')
        patient_description ['rawid'] = str
    except:
        print ("Could not parse the following ID as an AEH header: %s" % str)
        raise

    return patient_description, ('|'.join(bits[2:]) if len (bits) > 2  else None)

def parseRegExp (regexp):
    def parseHeader (str):
        try:
            bits = regexp.search(str.rstrip())
            patient_description = {}
            patient_description ['id']   =  bits.group(1)
            patient_description ['date']  = None
            patient_description ['rawid'] = str
        except:
            print ("Could not parse the following ID as the reg.exp. header: %s" % str)
            raise

        return patient_description, ('|'.join(bits[2:]) if len (bits.groups()) > 2  else None)
    return parseHeader

def parseLANL (str):
    try:
        bits = str.rstrip().split ('_')
        if len (bits) < 4:
            raise Exception ('Improperly formatted LANL header (need at least "subtype_country_accession_yyyy": %s' % str)

        patient_description = {}
        patient_description ['id']    =  bits[2]
        patient_description ['date']  = time.strptime (bits[3],'%Y')
        patient_description ['rawid'] = str
    except:
        print ("Could not parse the following ID as a LANL header: %s" % str)
        raise

    return patient_description, ('_'.join(bits[4:]) if len (bits) > 4  else None)

def parsePlain (str):

    patient_description = {}
    patient_description ['id']    = str
    patient_description ['date']  = None
    patient_description ['rawid'] = str

    return patient_description, None

def tm_to_datetime (tm_object):
    if tm_object is None: return None
    return datetime.datetime (tm_object.tm_year, tm_object.tm_mon, tm_object.tm_mday)

def datetime_to_tm (datetime_object):
    if datetime_object is None: return None
    return time.strptime(datetime_object.strftime("%Y-%m-%d"),"%Y-%m-%d")

def describe_vector (vector):
    vector.sort()
    l = len (vector)
    return {'count': l, 'min': vector[0], 'max': vector[-1], 'mean': sum(vector)/l, 'median':  vector [l//2] if l % 2 == 1 else 0.5*(vector[l//2-1]+vector[l//2]), "IQR": [vector [l//4], vector [(3*l)//4]] }
    
    
def _test_edge_support (triangles, sequence_file_name, hy_instance, p_value_cutoff):
    if hy_instance is None:
        hy_instance = hy.HyphyInterface ();
    script_path = os.path.realpath(__file__)
    hbl_path =  os.path.join(os.path.dirname(script_path), "data", "HBL", "TriangleSupport.bf")

    hy_instance.queuevar ('_py_sequence_file', sequence_file_name)
    triangle_spec = []
    for k in triangles:
        for i in k[:3]:
            triangle_spec.append (i)

    hy_instance.queuevar ('_py_triangle_sequences', triangle_spec)    
    hy_instance.runqueue (batchfile = hbl_path)
    if len(hy_instance.stderr):
        raise RuntimeError (hy_instance.stderr)


    return_object = [] # ((triangle), (p-values))

    for k, t in enumerate (triangles):
        return_object.append ( (t, hy_instance.getvar (str(k),hy.HyphyInterface.MATRIX)) )
    
    return return_object


#-------------------------------------------------------------------------------

class edge:
    def __init__ (self, patient1, patient2, date1, date2, visible, attribute = None, sequence_ids = None):
        if patient1 < patient2:
            self.p1    = patient1
            self.p2    = patient2
            self.date1   = date1
            self.date2   = date2
        else:
            self.p2    = patient1
            self.p1    = patient2
            self.date2   = date1
            self.date1   = date2

        if self.p1.id == self.p2.id:
            raise BaseException ("Trying to draw an edge between the same patient")
        self.visible = visible
        self.attribute = set()
        if attribute is not None:
            self.attribute.add (attribute)
        self.sequences = sequence_ids
        self.edge_reject_p = 0.
        self.is_unsupported = False
        
    def has_support (self):
        return self.is_unsupported == False

    def compute_direction (self, return_diff = False, min_days = 30, assume_missing_is_chronic = 180):
    # returns the node FROM which the edge is pointing AWAY
        if self.date1 and self.date2:
            if self.p2.edi:
                diff21 = (time.mktime(self.p2.edi) - time.mktime(self.date1))/(24*3600)
                if diff21 >= min_days:
                    return (self.p1, diff21) if return_diff else self.p1
                elif assume_missing_is_chronic is not None:
                    if self.p1.edi is None and diff21 >= -assume_missing_is_chronic:
                        return (self.p1, diff21) if return_diff else self.p1

            if self.p1.edi:
                diff12 = (time.mktime(self.p1.edi) - time.mktime(self.date2))/(24*3600)
                if diff12 >= min_days:
                    return (self.p2, diff12) if return_diff else self.p2
                elif assume_missing_is_chronic is not None:
                    if self.p2.edi is None and diff12 >= -assume_missing_is_chronic:
                        return (self.p2, diff12) if return_diff else self.p2

        return (None, 0) if return_diff else None

    def why_no_direction (self, min_days = 30):
        if self.date1 and self.date2:
            if self.p2.edi is None and self.p1.edi is None:
                return "No EDI"
            if self.p2.edi:
                diff21 = (time.mktime(self.p2.edi) - time.mktime(self.date1))/(24*3600)
                if diff21 < min_days:
                    if diff21 > 0:
                        return "Dates too close"
                    else:
                        return "Predates"
            if self.p1.edi:
                diff12 = (time.mktime(self.p1.edi) - time.mktime(self.date2))/(24*3600)
                if diff12 < min_days:
                    if diff12 > 0:
                        return "Dates too close"
                    else:
                        return "Predates"
        return "Missing dates"

    def direction(self, do_csv = False):
        dir = self.compute_direction()
        if dir and self.p1 == dir:
            return ["%s,%s,1"% (self.p1.id, self.p2.id)] if do_csv else ['"%s" -> "%s"' % (self.p1.id, self.p2.id), 'normal']
        elif dir and self.p2 == dir:
            return ["%s,%s,1"% (self.p2.id, self.p1.id)] if do_csv else ['"%s" -> "%s"' % (self.p2.id, self.p1.id), 'normal']

        return  ["%s,%s,0"% (self.p1.id, self.p2.id)] if do_csv else ['"%s" -> "%s"' % (self.p1.id, self.p2.id), 'none']

    def chrono_length_days (self):
        if self.date1 and self.date2:
            return abs(tm_to_datetime (self.date1) - tm_to_datetime (self.date2))
        return None

    def label (self):
        '''if self.date1 and self.date2:
            diff = self.chrono_length_days()
            return str (diff.days/7)'''
        return ''

    def __hash__ (self):
        return self.p1.__hash__() + self.p2.__hash__() + self.date1.__hash__() + self.date2.__hash__()

    def __comp__ (self,other):
        # 0: equal; 1: self is greater; -1: other is greater
        if self.p1 == other.p1 and self.p2 == other.p2:
            if self.date1 == other.date1 and self.date2 == other.date2:
                return 0
            if self.date1 is not None:
                if other.date1 is None:
                    return 1
                else:
                    if other.date1 > self.date1:
                        return -1
                    else:
                        if other.date1 < self.date1:
                            return 1
            else:
                if other.date1 is not None:
                    return -1

            if self.date2 is not None:
                if other.date2 is None:
                    return 1
                else:
                    if other.date2 > self.date2:
                        return -1
                    else:
                        if other.date2 < self.date2:
                            return 1
            return -1

        if self.p1 < other.p1:
            return -1
        elif self.p1 > other.p1:
            return 1
        elif self.p2 < other.p2:
            return -1
        return 1

    def update_attributes (self, desc):
        if desc is not None:
            self.attribute.add(desc)
            
    def update_sequence_info (self, seq_info):
        if seq_info is not None:
            self.sequences = seq_info

    def has_attribute (self, attr):
        return attr in self.attribute

    def remove_attribute (self, attr):
        self.attribute.discard (attr)

    def check_date (self, year, newer = False, weak = False):
        op = operator.__or__ if weak else operator.__and__
        if newer:
            return op(self.date1 == None or self.date1.tm_year >= year, self.date2 == None or self.date2.tm_year >= year)
        else:
            return op(self.date1 == None or self.date1.tm_year <= year ,self.date2 == None or self.date2.tm_year <= year)

    def check_exact_date (self, the_date, newer = False):
        if newer:
            return (self.date1 >= the_date) and (self.date2 >= the_date)
        else:
            return (self.date1 <= the_date) and (self.date2 <= the_date)

    def __lt__ (self, other):
        return self.__comp__ (other) == -1

    def __le__ (self, other):
        return self.__comp__ (other) <= 0

    def __gt__ (self, other):
        return self.__comp__ (other) == 1

    def __ge__ (self, other):
        return self.__comp__ (other) >= 0

    def __ne__ (self, other):
        return not self.__eq__ (other)

    def __eq__ (self,other):
        return self.p1 == other.p1 and self.p2 == other.p2 and self.date1 == other.date1 and self.date2 == other.date2

    def __repr__ (self):
        dir = self.compute_direction()
        if dir is None:
            dir = '--'
        elif dir == self.p1:
            dir = '->'
        else:
            dir = '<-'

        return "%s (%s) %s %s (%s)" % (self.p1.id, time.strftime ("%m-%d-%y",self.date1) if self.date1 is not None else 'None', dir, self.p2.id, time.strftime ("%m-%d-%y",self.date2) if self.date2 is not None else 'None')

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------

class patient:

    def __init__ (self, id):
        self.id                = id # a unique patient ID
        self.dates             = [] # date objects
        self.edi               = None # estimated date of infection
        self.stage             = 'Chronic' # disease stage
        self.treatment_date    = None # the date treatment started
        self.vl                = None # viral load at baseline
        self.degree            = 0
        self.cluster_id        = None
        self.naive             = None
        self.attributes        = set ()
        self.label             = None

    def __hash__ (self):
        return self.id.__hash__()

    def __comp__ (self,other):
        if self.id == other.id:
            return 0
        if self.id < other.id:
            return -1
        return 1

    def __eq__ (self,other):
        return self.id == other.id

    def __str__ (self):
        return "Patient %s (degree = %d, dates = %d, cluster_id = %s)" % (self.id, self.degree, len (self.dates), self.cluster_id)

    def __lt__ (self, other):
        return self.__comp__ (other) == -1

    def __le__ (self, other):
        return self.__comp__ (other) <= 0

    def __gt__ (self, other):
        return self.__comp__ (other) == 1

    def __ge__ (self, other):
        return self.__comp__ (other) >= 0

    def __ne__ (self, other):
        return not self.__eq__ (other)

    def __repr__ (self):
        return self.__str__()

    def add_attribute (self, attrib):
        if attrib is not None:
            self.attributes.add (attrib)

    def remove_attribute (self, attrib):
        self.attributes.discard (attrib)

    def has_attribute (self, attr):
        return attr in self.attributes

    def add_date (self, date):
        if date not in self.dates:
            self.dates.append (date)

    def add_degree (self):
        self.degree += 1

    def add_edi  (self, edi):
        self.edi = edi

    def get_edi (self):
        return self.edi

    def add_stage  (self, stage):
        self.stage = stage

    def add_treatment  (self, drugz):
        self.treatment_date = drugz

    def get_label (self):
        return self.label

    def set_label (self, l):
        self.label = l

    def add_vl (self, vl, date = None):
        if self.vl is None:
            self.vl = []

        for k in range (len (self.vl)):
            if self.vl[k][0] is not None and self.vl[k][0] > date:
                self.vl.insert (k, [date, float(vl)])
                return None


        self.vl.append ([date, float(vl)])


    def get_vl_by_date (self, date):
        if self.vl is not None:
            if len  (self.vl):
                idx = bisect_left ([d[0] for d in self.vl], date)
                if idx == len (self.vl):
                    return self.vl [idx-1]
                if idx > 0 and idx < len (self.vl):
                    if (date-self.vl[idx-1][0] < self.vl[idx][0]-date):
                        return self.vl [idx-1]

                return self.vl [idx]

        return None

    def get_vl (self, date = None):
        if self.vl is not None and len (self.vl):
            return sum ([v[1] for v in self.vl]) / float (len (self.vl))

        return None

    def add_naive (self, naive):
        self.naive = naive

    def get_followup_length (self, date):
        if self.dates[0] is None:
            return None

        return date - tm_to_datetime (min (self.dates))

    def get_baseline_date (self, complete = False):
        if self.dates[0] is None:
            return None

        if complete:
            return min (self.dates)
        return min ([k.tm_year for k in self.dates if k is not None])

    def get_latest_date (self, complete = False):
        if complete:
            return max (self.dates)
        return max ([k.tm_year for k in self.dates if k is not None])

    def get_sample_count (self):
        return len(self.dates)

    def get_length_of_followup (self):
        if None not in self.dates:
            d1 = tm_to_datetime (self.dates[0])
            if len (self.dates) > 1:
                self.dates.sort()
                d2 = tm_to_datetime(self.dates[-1])
                return d2 - d1
        return datetime.timedelta (0)


    def get_treatment_since_edi (self):
        if self.treatment_date != None and self.edi != None and self.treatment_date >= self.edi:
            d1 = tm_to_datetime (self.treatment_date)
            d2 = tm_to_datetime (self.edi)
            return d1-d2
        return None

    def get_dot_string (self, year_vis = None):

        '''if self.id [0:4] == '0501':
            lab = self.id [4:]
        else:
            lab = self.id

        return '"%s" [label = "%s"];\n' % (self.id, lab) '''

        shape = 'ellipse'
        color = 'white'
        label = str(get_vl()) if self.get_vl() is not None else ""

        edi_info = self.get_treatment_since_edi()


        if edi_info:
            if edi_info.days <= 30:
                color = 'green'
            else:
                color = 'yellow'
            #label = str(edi_info.days/7)

        if self.naive:
            color = 'red'

        if year_vis is not None:
            if self.get_baseline_date () > year_vis:
                return '"%s" [label = "%s", fillcolor = "%s", shape = %s, style = "invis"];\n' % (self.id, label , color, shape)

        return '"%s" [label = "%s", fillcolor = "%s", shape = %s];\n' % (self.id, label , color, shape)

#-------------------------------------------------------------------------------

class transmission_network:

    def __init__ (self):
        self.nodes     = {}
        self.edges     = {}
        self.adjacency_list = None
        self.sequence_ids = {} # this will store unique sequence ids keyed by edge information (pid and date)

    def read_from_csv_file (self,file_name, formatter = None, distance_cut = None, default_attribute = None, bootstrap_mode = False):
        if formatter is None:
            formatter = parseAEH
        edgeReader = csv.reader(file_name)
        header = next(edgeReader)
        edgeAnnotations = {}
        if len (header) < 3:
            raise IOError ('transmission_network.read_from_csv_file() : Expected a .csv file with at least 3 columns as input')

        handled_ids = set ()

        for line in edgeReader:
            distance = float(line[2])
            if distance_cut is not None and distance > distance_cut:
                self.ensure_node_is_added (line[0], formatter, default_attribute, bootstrap_mode, handled_ids)
                self.ensure_node_is_added (line[1], formatter, default_attribute, bootstrap_mode, handled_ids)
                continue
            edge = self.add_an_edge(line[0],line[1],distance,formatter,default_attribute,bootstrap_mode)
            if edge is not None and len (line) > 3:
                 edgeAnnotations[edge] = line[2:]

        return edgeAnnotations

    def ensure_node_is_added (self, id1, header_parser, default_attribute, bootstrap_mode, cache):
        if id1 not in cache:
            cache.add (id1)
            if header_parser == None:
                header_parser = parseAEH

            patient1,attrib = header_parser (id1)
            self.insert_patient (patient1['id'],patient1['date'], False, attrib)

    def sample_from_network (self, how_many_nodes = 100, how_many_edges = None, node_sampling_bias = 0.0):
        if how_many_edges is not None:
            if how_many_edges >= len (self.edges):
                return self
            subset_network = transmission_network()
            sampled_edges = random.sample (list(self.edges), how_many_edges)
            for an_edge in sampled_edges:
                subset_network.add_an_edge (an_edge.p1.id, an_edge.p2.id, self.edges[an_edge], header_parser = parsePlain)
            return subset_network


        if how_many_nodes >= len (self.nodes):
            return self

        subset_network = transmission_network()

        if node_sampling_bias > 0.:
            nodes     = set ()
            left_over = set (list(self.nodes))
            connected_to_existing = set ()

            node_neighbs   = {}
            sampling_probs = {}

            for a_node in self.nodes:
                node_neighbs [a_node]   = self.get_node_neighborhood (a_node.id,True,False)
                sampling_probs [a_node] = 1.

            nodes.add                    (a_node)
            left_over.remove             (a_node)
            connected_to_existing.update (node_neighbs[a_node])
            connected_to_existing.add (a_node)

            for a_node in connected_to_existing:
                sampling_probs[a_node] += node_sampling_bias

            upper_bound   = len(left_over) + (len(connected_to_existing)-1)*node_sampling_bias

            while len (nodes) < how_many_nodes:
                up_to = random.random() * upper_bound
                local_sum = 0
                for a_node in left_over:
                    local_sum += sampling_probs[a_node]
                    if local_sum >= up_to:
                        break

                if local_sum < up_to:
                    raise Exception ('Sampling fubar')

                nodes.add                    (a_node)
                connected_to_existing.add    (a_node)
                new_connections = node_neighbs[a_node]-connected_to_existing
                connected_to_existing.update (new_connections)
                #print (a_node.id, len(new_connections), upper_bound, up_to, sum ([sampling_probs[k] for k in left_over]))
                left_over.remove             (a_node)

                upper_bound   += len(new_connections)*node_sampling_bias - sampling_probs[a_node]

                for a_node in new_connections:
                    sampling_probs[a_node] += node_sampling_bias



            #print (sampling_probs.values())
        else:
            nodes = random.sample (list(self.nodes), how_many_nodes)

        for a_node in nodes:
            this_node = self.nodes[a_node]
            subset_network.insert_patient (this_node.id, this_node.dates[0], False, None)

        for an_edge in self.edges:
            if an_edge.p1 in subset_network.nodes and an_edge.p2 in subset_network.nodes:
                subset_network.add_an_edge (an_edge.p1.id, an_edge.p2.id, self.edges[an_edge], header_parser = parsePlain)

        return subset_network

    def create_a_random_network (self, network_size = 100):
        self.insert_patient (1, None, False, None)
        for node_id in range (2, network_size+1):
            self.add_an_edge (node_id, random.randint (1, node_id-1), 1, header_parser = parsePlain)

        return self

    def add_contemporaneuos_edges (self, diff, prob):
        self.compute_clusters()
        clusters = self.retrieve_clusters()
        for cluster_id in clusters.keys():
            if cluster_id is not None:
                pairs = itertools.combinations(clusters[cluster_id],2)
                for node_pair in pairs:
                     d1 = tm_to_datetime (node_pair[0].dates[0])
                     d2 = tm_to_datetime (node_pair[1].dates[0])
                     #print (d1, d2, diff)
                     if (d1 < d2):
                        if d2-d1 >= diff: continue
                     else:
                        if d1-d2 >= diff: continue

                     if random.random () < prob:
                        #print (node_pair[0].dates[0], node_pair[1].dates[0])
                        self.add_an_edge ("|".join ([node_pair[0].id,time.strftime("%m%d%Y",node_pair[0].dates[0])]), "|".join ([node_pair[1].id,time.strftime("%m%d%Y",node_pair[1].dates[0])]), 0.01, header_parser = parseAEH)

    def create_a_pref_attachment_network (self, network_size = 100, random_attachment = 0.0, start_new_tree = 0.0, start_date = None, tick_rate = None):
        attach_to = [1,2]
        current_date = start_date

        self.insert_patient (str(1), datetime_to_tm(current_date), False, None)
        self.insert_patient (str(2), datetime_to_tm(current_date), False, None)

        for node_id in range (3, network_size+1):
            if current_date is not None:
                current_date += datetime.timedelta (days = random.expovariate (tick_rate))
                #print (current_date)
            if start_new_tree > 0.0 and random.random () < start_new_tree:
                self.insert_patient (str(node_id), datetime_to_tm(current_date), False, None)
                attach_to.append (node_id)
                continue


            if random_attachment > 0.0 and random.random() < random_attachment:
                k = random.randint (1,node_id-1)
            else:
                k = random.choice (attach_to)

            if current_date is not None:
                self.add_an_edge ("|".join([str(node_id),current_date.strftime("%m%d%Y")]), "|".join ([str(k),time.strftime("%m%d%Y",self.has_node_with_id(str(k)).dates[0])]), 1, header_parser = parseAEH)
            else:
                self.add_an_edge (str(node_id), str(k), 1, header_parser = parsePlain)
            attach_to.extend ([k, node_id])

        return self

    def insert_patient (self, id, date, add_degree, attributes):
        pat = patient (id)
        if pat not in self.nodes:
            self.nodes[pat] = pat

        pat = self.nodes[pat]
        pat.add_date  (date)
        if add_degree:
            pat.add_degree()

        pat.add_attribute (attributes)
        return pat

    def make_sequence_key (self, id, date):
        if date != None:
            return "-".join((id,time.strftime ("%m-%d-%Y",date)))
        return id

    def add_edi  (self, edi):
        for node in self.nodes:
            if node.id in edi:
                #[geno_date, drug_date, edi_date, viral_load, naive]
                node.add_treatment (edi[node.id][1])
                node.add_stage (edi[node.id][2])
                node.add_edi (edi[node.id][3])
                if edi[node.id][3] is not None:
                    node.add_stage (edi[node.id][2])

                node.add_vl (edi[node.id][4])
                node.add_naive (edi[node.id][5])

    def add_edi_json  (self, edi):
        for pid in edi:
            p = patient (pid)
            if p in self.nodes:
                node = self.nodes[p]
                edi_record = edi[pid]
                if 'Stage' in edi_record:
                    node.add_stage (edi_record['Stage'])
                if 'EDI' in edi_record:
                    node.add_edi (edi_record['EDI'])
                if 'VL' in edi_record:
                    for vl_record in edi_record['VL']:
                        node.add_vl (vl_record[1], vl_record[0])
                if 'ARV' in edi_record:
                    node.add_treatment (edi_record['ARV'])


    def clustering_coefficients (self, node_list = None):
        clustering_coefficiencts = {}

        if self.adjacency_list is None:
            self.compute_adjacency ()

        if node_list is None:
            node_list = self.adjacency_list.keys()

        precompute_neighborhoods = {}
        for a_node in node_list:
            precompute_neighborhoods[a_node] = set(self.get_node_neighborhood (a_node.id))

        for a_node in node_list:
            if a_node in self.adjacency_list:
                my_hbhd = precompute_neighborhoods[a_node]
                nbhd_size = len (my_hbhd)
                if nbhd_size > 1:
                    edges_found = 0
                    for nb_node in my_hbhd:
                        edges_found += len(precompute_neighborhoods[nb_node].intersection (my_hbhd))
                    clustering_coefficiencts [a_node] = float(edges_found) / nbhd_size / (nbhd_size-1)
        return clustering_coefficiencts

    def randomize_attribute (self, attribute_value, clusters = None):
        if clusters is None:
            partition = [self.nodes,]
        else:
            partition = clusters

        for subset in partition:
            count = 0
            for node in subset:
                if node.has_attribute (attribute_value):
                    count += 1
                    node.remove_attribute (attribute_value)

            for node in random.sample (list(subset), count):
                node.add_attribute (attribute_value)

    def edges_sharing_an_attribute (self, attribute_value = None, reduce_edges = True, ignore_visible = False):
        #if attribute_value == None, then any shared attributes count

        result = {'compared' : 0, 'shared' : 0};

        for anEdge in self.edges if reduce_edges == False else self.reduce_edge_set():
            if anEdge.visible or ignore_visible:
                result['compared'] += 1
                if attribute_value is None:
                    result['shared'] += (1 if (len(anEdge.p1.attributes.intersection (anEdge.p2.attributes)) > 0) else 0)
                else:
                    result['shared'] += (1 if anEdge.p1.has_attribute (attribute_value) and anEdge.p2.has_attribute (attribute_value) else 0)

        return result;

    def has_an_edge        (self, id1, id2):
        test_edge = edge (patient (id1), patient (id2), None, None, True)
        return self.edges[test_edge] if test_edge in self.edges else None

    def has_node_with_id        (self, id):
        pat_with_id = patient (id)
        if pat_with_id in self.nodes:
            return self.nodes[pat_with_id]
        return None

    def get_all_edges_linking_to_a_node  (self, id1, ignore_visible = False, use_direction = False, incoming = False, add_undirected = False, only_undirected = False, reduce_edges = True):
        list_of_nodes = set()
        pat = patient (id1)
        for anEdge in self.edges if reduce_edges == False else self.reduce_edge_set():
            if anEdge.visible or ignore_visible:
                if pat == anEdge.p2 or pat == anEdge.p1:
                    if use_direction:
                        dir = anEdge.compute_direction()
                        if dir is not None:
                            if only_undirected:
                                continue
                            if ((not incoming and dir != pat) or (incoming and dir == pat)):
                                continue
                        elif not add_undirected and not only_undirected:
                            continue

                    list_of_nodes.add (anEdge)
        return list_of_nodes

    def get_node_neighborhood  (self, id1, ignore_visible = False, use_direction = False, incoming = False, add_undirected = False, only_undirected = False):
        list_of_nodes = set()
        list_of_edges = self.get_all_edges_linking_to_a_node(id1,ignore_visible,use_direction, incoming, add_undirected, only_undirected)
        pat = patient (id1)
        for anEdge in list_of_edges:
            if anEdge.p1 == pat:
                list_of_nodes.add (anEdge.p2)
            else:
                list_of_nodes.add (anEdge.p1)

        return list_of_nodes

    def summarize_bootstrap (self):
        edge_support = {}
        for edge in self.edges:
            pair = (edge.p1, edge.p2)
            if pair not in edge_support:
                edge_support[pair] = set()
            edge_support[pair].update (edge.attribute)

        for k in edge_support:
            edge_support[k] = len (edge_support[k])

        return edge_support


    def add_an_edge (self, id1, id2, distance, header_parser = None, edge_attribute = None, bootstrap_mode = False, node_only = False):
        if header_parser == None:
            header_parser = parseAEH

        patient1,attrib = header_parser (id1)
        patient2,attrib = header_parser (id2)
        same = patient1['id'] == patient2['id']

        p1 = self.insert_patient (patient1['id'],patient1['date'], not same and not node_only, attrib)
        p2 = self.insert_patient (patient2['id'],patient2['date'], not same and not node_only, attrib)

        pid1 = self.make_sequence_key (patient1['id'],patient1['date'])
        if pid1 not in self.sequence_ids:
            self.sequence_ids [pid1] = patient1["rawid"]
            
        pid2 = self.make_sequence_key (patient2['id'],patient2['date'])
        if pid2 not in self.sequence_ids:
            self.sequence_ids [pid2] = patient2["rawid"]

        if node_only == False:
            if not same:
                new_edge = edge (p1,p2,patient1['date'],patient2['date'],True, edge_attribute, (patient1["rawid"], patient2["rawid"]))
                if new_edge not in self.edges:
                    if not bootstrap_mode or edge_attribute is None:
                        self.edges [new_edge] = distance
                else:
                    for k in self.edges:
                        if k == new_edge:
                            k.update_attributes (edge_attribute)
                            k.update_sequence_info ((patient1["rawid"], patient2["rawid"]))
                            break

                    if distance < self.edges [new_edge]:
                        self.edges [new_edge] = distance

                return new_edge

        return None

    def compute_adjacency (self,edges=False,edge_set = None, both = False):
        self.adjacency_list = {}
        
        for anEdge in (edge_set if edge_set is not None else self.edges):
            if anEdge.visible:
                if anEdge.p1 not in self.adjacency_list: self.adjacency_list [anEdge.p1] = set()
                if anEdge.p2 not in self.adjacency_list: self.adjacency_list [anEdge.p2] = set()
                if edges:
                    # check for duplication
                    processed = False
                    for an_edge in self.adjacency_list [anEdge.p1]:
                        if an_edge.p1 == anEdge.p1 and an_edge.p2 == anEdge.p2:
                            if an_edge > anEdge: #existing is "greater", replace
                                self.adjacency_list[anEdge.p1].remove (an_edge)
                                self.adjacency_list[anEdge.p2].remove (an_edge)
                            else:
                                processed = True
                            break

                    if not processed:
                        self.adjacency_list [anEdge.p1].add (anEdge)
                        self.adjacency_list [anEdge.p2].add (anEdge)
                elif both:
                    # check for duplication
                    processed = False
                    for a_node, an_edge in self.adjacency_list [anEdge.p1]:
                        if an_edge.p1 == anEdge.p1 and an_edge.p2 == anEdge.p2:
                            if an_edge > anEdge: #existing is "greater", replace
                                self.adjacency_list[anEdge.p1].remove ((a_node,an_edge))
                                self.adjacency_list[anEdge.p2].remove ((a_node,an_edge))
                            else:
                                processed = True
                            break

                    if not processed:
                        self.adjacency_list [anEdge.p1].add ((anEdge.p2,anEdge))
                        self.adjacency_list [anEdge.p2].add ((anEdge.p1,anEdge))
                
                else:
                    self.adjacency_list [anEdge.p1].add (anEdge.p2)
                    self.adjacency_list [anEdge.p2].add (anEdge.p1)


    def compute_path_stat (self, distances, stat = "mean"):
        result = {}
        node_count = len (distances['ordering'])
        for i in range(node_count):
            d = 0
            for j,p in enumerate(distances['distances'][i]):
                if p is None:
                    if j==i: continue
                    d = None
                    break
                else:
                    d += p

            result[distances['ordering'][i]] = d / (node_count-1)

        return result


    def compute_shortest_paths (self, subset = None, use_actual_distances = False):
        self.compute_adjacency()

        if subset is None:
            subset = self.adjacency_list.keys()

        node_count = len (subset)
        distances  = []

        for a_node in (subset):
            distances.append ([None for k in range (node_count)])

        for index,a_node in enumerate(subset):
            for index2, second_node in enumerate(subset):
                if second_node != a_node:
                    if second_node in self.adjacency_list [a_node]:
                        distances[index][index2] = 1
                        distances[index2][index] = 1


        distances2 = deepcopy (distances)

        for index_k, n_k in enumerate (subset):
            for index_i,n_i in enumerate(subset):
                for index_j, n_j in enumerate(subset):
                    if n_i != n_j:
                        d_ik = distances[index_k][index_i]
                        d_jk = distances[index_k][index_j]
                        d_ij = distances[index_i][index_j]
                        if d_ik is not None and d_jk is not None:
                            d_ik += d_jk
                            if d_ij is None or d_ij > d_ik:
                                distances2 [index_i][index_j] = d_ik
                                distances2 [index_j][index_i] = d_ik
                                continue
                        distances2 [index_j][index_i] = distances [index_j][index_i]
                        distances2 [index_i][index_j] = distances [index_i][index_j]

            t = distances2
            distances2 = distances
            distances = t

        return {'ordering': subset, 'distances': distances}

    def get_all_treated_within_range (self, daterange, outside = False):
        selection = []
        for node in self.nodes:
            tedi = node.get_treatment_since_edi()
            if tedi and (tedi > daterange if outside else tedi <= daterange):
                selection.append (node)
        return selection

    def get_all_naive (self):
        selection = []
        for node in self.nodes:
            if node.naive:
                selection.append (node)
        return selection

    def report_multiple_samples (self, minfo):
        counts = describe_vector ([k[0] for k in minfo])
        fup    = describe_vector ([k[1] for k in minfo])

        return {'count': counts['count'], 'samples' : counts, 'followup' : fup}


    def get_edge_node_count (self, attributes_to_check = None):
        vis_nodes = set ()
        edge_set  = set ()
        multiple_samples = []
        edge_count = 0

        nodes_by_stage = {}

        for edge in self.edges:
            if edge.visible:
                edge_count += 1
                for p in [edge.p1,edge.p2]:
                    if p not in vis_nodes:
                        if attributes_to_check is not None:
                            if not attributes_to_check.issubset (p.attributes):
                                continue
                        vis_nodes.add (p)
                        if p.stage not in nodes_by_stage:
                            nodes_by_stage [p.stage] = 1
                        else:
                            nodes_by_stage [p.stage] += 1

                        if p.get_sample_count() > 1:
                            multiple_samples.append([p.get_sample_count(),p.get_length_of_followup()])


                edge_set.add ((edge.p1, edge.p2))


        return {'edges': len(edge_set), 'nodes': len(vis_nodes), 'total_edges': edge_count, 'multiple_dates':[[k[0],k[1].days] for k in multiple_samples], 'total_sequences': len(vis_nodes) + sum([k[0] for k in multiple_samples]) - len(multiple_samples), 'stages' : nodes_by_stage }

    def clear_adjacency (self, clear_filter = True):
        if self.adjacency_list is not None:
            del self.adjacency_list
            self.adjacency_list = None
            if clear_filter:
                self.clear_filters()


    def apply_disease_stage_filter  (self, stages, do_clear = True, do_exclude = False):
        if do_clear : self.clear_adjacency()
        vis_count = 0
        for edge in self.edges:
            if edge.visible:
                if do_exclude:
                    edge.visible = edge.p1.stage not in stages and edge.p2.stage not in stages
                else:
                    edge.visible = edge.p1.stage in stages and edge.p2.stage in stages
                vis_count += edge.visible

        return vis_count

    def apply_date_filter     (self, edge_year, newer = False, do_clear = True):
        if do_clear : self.clear_adjacency()
        vis_count = 0
        for edge in self.edges:
            if edge.visible:
                edge.visible = edge.check_date (edge_year, newer)
                vis_count += edge.visible
        return vis_count

    def apply_exact_date_filter     (self, the_date, newer = False, do_clear = True):
        if do_clear : self.clear_adjacency()
        vis_count = 0
        for edge in self.edges:
            if edge.visible:
                edge.visible = edge.check_exact_date (the_date, newer)
                vis_count += edge.visible
        return vis_count

    def apply_distance_filter (self, distance, do_clear = True):
        if do_clear : self.clear_adjacency()
        vis_count = 0
        for edge in self.edges:
            if edge.visible:
                edge.visible = self.edges[edge] <= distance
                vis_count += edge.visible
        return vis_count

    def apply_id_filter (self, list, strict = False, do_clear = True):
        if do_clear : self.clear_adjacency()
        vis_count = 0
        for edge in self.edges:
            if edge.visible:
                if strict:
                    edge.visible = edge.p1.id in list and edge.p2.id in list
                else:
                    edge.visible = edge.p1.id in list or edge.p2.id in list
                vis_count += edge.visible
        return vis_count

    def apply_removed_edge_filter (self, do_clear = True):
        if do_clear : self.clear_adjacency()
        vis_count = 0
        for edge in self.edges:
            if edge.visible:
                edge.visible = edge.has_support()
            vis_count += edge.visible
        return vis_count

    def apply_attribute_filter (self, attribute_value, do_clear = True, strict = False):
        if do_clear : self.clear_adjacency()
        vis_count = 0
        for edge in self.edges:
            if edge.visible:
                if strict:
                    edge.visible = edge.p1.has_attribute(attribute_value) and edge.p2.has_attribute(attribute_value)
                else:
                    edge.visible = edge.p1.has_attribute(attribute_value) or edge.p2.has_attribute(attribute_value)
                vis_count += edge.visible
        return vis_count

    def apply_cluster_filter (self, cluster_ids, exclude = True, do_clear = True): # exclude all sequences in a given cluster(s)
        if do_clear : self.clear_adjacency()
        vis_count = 0

        for edge in self.edges:
            if edge.visible:
                if edge.p1.cluster_id in cluster_ids or edge.p2.cluster_id in cluster_ids:
                    edge.visible = not exclude
                else:
                    edge.visible = exclude

            vis_count += edge.visible

        return vis_count


    def retrieve_clusters (self, singletons = True):
        clusters = {}
        for node in self.nodes:
            #if node.cluster_id == None:
            #    raise BaseException ('Called return_clusters but node %s had no associated cluster ID' % node)
            if node.cluster_id not in clusters:
                clusters [node.cluster_id] = []
            clusters [node.cluster_id].append (node)

        if not singletons:
            clusters.pop (None,None)
        return clusters


    def clear_filters           (self):
        for edge in self.edges:
            edge.visible = True

    def cluster_size_by_node (self):
        if self.adjacency_list == None:
            self.compute_adjacency ()
        self.compute_clusters()
        clusters     = self.retrieve_clusters ()
        size_by_node = {}
        for c in clusters:
            if c is not None:
                for node in clusters[c]:
                    size_by_node[node] = len(clusters[c])
        for node in self.nodes:
            if node not in size_by_node:
                size_by_node[node] = 1
        return size_by_node


    def drop_singleton_nodes (self):
        if self.adjacency_list == None:
            self.compute_adjacency ()

        drop_these = set ()
        for a_node in self.nodes:
            if a_node not in self.adjacency_list:
                drop_these.add (a_node)

        for delete_me in drop_these:
            del self.nodes[delete_me]

    def compute_clusters (self, singletons = False):
        if self.adjacency_list == None:
            self.compute_adjacency ()

        for aNode in self.nodes:
            aNode.cluster_id = None

        cluster_id = [0] # this will pass the object by reference

        for node in self.nodes:
            if (singletons or node in self.adjacency_list) and node.cluster_id == None:
                self.breadth_first_traverse (node, cluster_id)

    def breadth_first_traverse (self, node, cluster_id):
        if node.cluster_id == None:
            cluster_id [0] += 1
            node.cluster_id = cluster_id [0]
        if node in self.adjacency_list:
            for neighbor_node in self.adjacency_list[node]:
                if neighbor_node.cluster_id == None:
                    neighbor_node.cluster_id = node.cluster_id
                    self.breadth_first_traverse(neighbor_node, cluster_id)

    def generate_csv (self, file):
        file.write ("ID1,ID2,Distance")
        for edge in self.edges:
            if edge.visible:
                file.write("%s,%s,%g\n" % (edge.p1.id, edge.p2.id, self.edges[edge]))

    def write_clusters (self, file):
        file.write ("SequenceID,ClusterID\n")
        for node in self.nodes:
            if node.cluster_id != None:
                file.write ("%s,%d\n" % (self.sequence_ids[self.make_sequence_key (node.id, node.dates[0])],node.cluster_id))


    def reduce_edge_set (self, attribute_merge = True):
        byPairs = {}
        for anEdge in self.edges:
            patient_pair = (anEdge.p1, anEdge.p2)
            if patient_pair in byPairs:
                byPairs[patient_pair].append (anEdge)
            else:
                byPairs[patient_pair] = [anEdge]

        edge_set = set ()
        for patient_pair in byPairs:
            representative_edge = min (byPairs[patient_pair])
            if attribute_merge:
                attribute_set = set ()
                for an_edge in byPairs[patient_pair]:
                    attribute_set = attribute_set.union (an_edge.attribute)
                for attr in attribute_set:
                    representative_edge.update_attributes(attr)

            edge_set.add (representative_edge)
        return edge_set

    def prune_all_edges_lacking_support (self):
        byPairs = {}
        for anEdge in self.edges:
            patient_pair = (anEdge.p1, anEdge.p2)
            if patient_pair in byPairs:
                byPairs[patient_pair].append (anEdge)
            else:
                byPairs[patient_pair] = [anEdge]

        for patient_pair in byPairs:
            representative_edge = min (byPairs[patient_pair])
            if representative_edge.has_support () == False:
                for e in byPairs[patient_pair]:
                    del self.edges[e]
                
    

    def generate_dot (self, file, year_vis = None, reduce_edges = True):

        if self.adjacency_list is None:
             self.compute_adjacency()

        file.write ('digraph G { overlap="voronoi";\n outputorder = edgesfirst;\nnode[style=filled];\n');
        nodes_drawn = {}

        directed = {'undirected':0, 'directed':0}

        for edge in self.edges if reduce_edges == False else self.reduce_edge_set():
            if edge.visible:
                distance = self.edges[edge]

                if edge.p1 not in nodes_drawn:
                    nodes_drawn[edge.p1] = edge.p1.get_baseline_date()
                    file.write (edge.p1.get_dot_string(year_vis))
                if edge.p2 not in nodes_drawn:
                    nodes_drawn[edge.p2] = edge.p2.get_baseline_date()
                    file.write (edge.p2.get_dot_string(year_vis))

                if isinstance(edge.compute_direction(),type(None)):
                    directed ['undirected'] += 1
                else:
                    directed ['directed'] += 1
                edge_attr = edge.direction()

                if year_vis is not None:
                    if edge.check_date (year_vis) == False:
                        file.write ('%s [style="invis" arrowhead = "%s"];\n' % (edge_attr[0], edge_attr[1]));
                        continue


                file.write ('%s [style="bold" label = "%s" arrowhead = "%s"];\n' % (edge_attr[0], edge.label(), edge_attr[1]));

        file.write ("\n};")
        return directed

    def generate_delimited (self, file, year_vis = None, reduce_edges = True):

        if self.adjacency_list is None:
             self.compute_adjacency()

        file.write("%s\n"%','.join(['ID1','ID2','Linktype']))
        for edge in self.edges if reduce_edges == False else self.reduce_edge_set():
            if edge.visible:
                distance = self.edges[edge]


                edge_attr = edge.direction(do_csv=True)

                if year_vis is not None:
                    if edge.check_date (year_vis) == False:
                        continue

                file.write ('%s\n' % (edge_attr[0]));


    def spool_pairwise_distances (self,file,baseline = False):
        file.write (','.join (['Seq1','Seq2','Distance']))
        file.write ('\n')
        for ext_edge in self.edges:
            if baseline:
                 if ext_edge.p1.get_baseline_date (True) != ext_edge.date1 or ext_edge.p2.get_baseline_date (True) != ext_edge.date2:
                    continue
            file.write (','.join ([ext_edge.p1.id, ext_edge.p2.id, str(self.edges[ext_edge])]))
            file.write ('\n')

    def get_node_degree_list (self, year_cap = None, do_direction = False, id_list = None, attribute_selector = None, clear_filters = True):
        degree_list = {}
        self.clear_adjacency(clear_filter = clear_filters)
        self.compute_adjacency(do_direction)

        if id_list:
            id_list = [self.has_node_with_id(k) for k in id_list]
        else:
            if attribute_selector is not None:
                id_list = [k for k in self.nodes if k.has_attribute(attribute_selector)]
            else:
                id_list = self.nodes

        for node in id_list:
            if year_cap is not None and node.get_baseline_date() > year_cap:
                    degree_list[node] = None
            else:
                if node in self.adjacency_list:
                    if do_direction:
                        degs = [0,0,0,0] # undir, out-edges, in-edges
                        for e in self.adjacency_list[node]:
                            dir = e.compute_direction()
                            if dir is None:
                                degs[0] += 1
                            elif dir == node:
                                degs[1] += 1
                            else:
                                degs[2] += 1
                        degs[3] = sum (degs[:3])
                        degree_list[node] = degs
                    else:
                        degree_list[node] = len (self.adjacency_list[node])
                else:
                    degree_list[node] = 0 if not do_direction else [0,0,0,0]

        return degree_list

    def sample_subset    (self, size, filter_attribute = None, use_connected_nodes = False):
        if use_connected_nodes:
            self.compute_adjacency ()
            if filter_attribute is not None:
                random.sample ([n for n in self.adjacency_list if n.has_attribute(filter_attribute)], int (size))
            return random.sample (list(self.adjacency_list), int (size))
        if filter_attribute is not None:
            random.sample ([n for n in self.nodes if n.has_attribute(filter_attribute)], int (size))
        return random.sample (list(self.nodes), int (size))


    def generate_random_edges (self, edge_count, only_new = True, node_set = None, default_attr = None, distance = 0.01, use_preferential_attachment = False):
        edges_added = set ()
        added = 0
        if node_set is None:
            node_set = list(self.nodes)

        if use_preferential_attachment:
            sample_from = list ()
            if self.adjacency_list is None:
                self.compute_adjacency ()
            for a_node in node_set:
                upto = max(1,len(self.adjacency_list[a_node])) if a_node in self.adjacency_list else 1
                for k in range (upto):
                    sample_from.append (a_node)
        else:
            sample_from = node_set


        while added < edge_count:
            n1, n2 = random.sample (sample_from, 2)
            if n1 == n2: continue
            added_edge = self.add_an_edge (n1.id, n2.id, distance, header_parser = parsePlain, edge_attribute = default_attr)
            added += 1 if (added_edge is not None or not only_new) else 0
            if added_edge:
                edges_added.add (added_edge)

        return edges_added

    def delete_edge_subset (self, edges):
        for an_edge in edges:
            if an_edge in self.edges:
                del self.edges[an_edge]

    def sample_subset_year_list    (self, years):
        selected_nodes = set()
        if len (years) > len (self.nodes):
            return None
        for a_year in years:
            a_sample = random.sample (list(self.nodes), 1)[0]
            while (a_sample.get_baseline_date() != a_year or a_sample in selected_nodes):
                a_sample = random.sample (list(self.nodes), 1)[0]
            selected_nodes.add (a_sample)
        return selected_nodes


    def output_sequence_names (self):
        pass

    def type_of_adjacency_list (self):
        if self.adjacency_list and self.adjacency_list is not None:
            if isinstance ([(k, i) for k, i in enumerate (self.adjacency_list) if k == 0][0][1], patient):
                return 'patient'
            else:
                return 'edge'
        return None

    def find_all_triangles (self, edge_set):
        triangles                   =  set ()
        #sequences_involved_in_links =  set ()
        #sequence_pairs              =  set ()


        if self.adjacency_list == None:
            self.compute_adjacency (both = True, edge_set = edge_set)
            
        print ("Locating triangles in the network", file = sys.stderr)

        '''
        for an_edge in edge_set:
            if an_edge.sequences is not None:
                sequence_pairs.add (an_edge.sequences)
                for seq in an_edge.sequences:
                    sequences_involved_in_links.add (seq)
            
        for this_pair in sequence_pairs:
            for third_apex in sequences_involved_in_links: 
                if (((third_apex, this_pair[0]) in sequence_pairs or (this_pair[0],third_apex) in sequence_pairs)
                    and ((third_apex, this_pair[1]) in sequence_pairs or (this_pair[1],third_apex) in sequence_pairs)):
                    triangles.add ((this_pair[0], this_pair[1], third_apex))
        '''
        
        # create a list of the form 
        # node[id] = list of 
        #   [node_id] = edge_id
        
        adjacency_map = {}
        for node, edge_list in self.adjacency_list.items():
            node_neighborhood = {}
            for n, e in edge_list:
                node_neighborhood [n] = e
            adjacency_map[node] = node_neighborhood
        
        triangle_nodes = set ()
        triangle_nodes_all = set ()
        
        count_by_sequence = {}
        
        for node, neighbors in adjacency_map.items():
            if len (neighbors) > 1: # something to do
                for node2 in neighbors:
                    for node3 in adjacency_map[node2]: 
                        if node in adjacency_map[node3]:   
                            triad = sorted([node, node2, node3])
                            triad = (triad[0], triad[1], triad[2])
                            if triad not in triangle_nodes:
                                sequence_set = set ()
                                for triangle_edge in [adjacency_map[node][node2], adjacency_map[node][node3], adjacency_map[node2][node3]]:
                                    sequence_set.update (triangle_edge.sequences)
                                if  len (sequence_set) == 3:
                                    triangle_nodes.add (triad)
                                    seqs = [k for k in sequence_set]
                                    triangles.add ((seqs[0], seqs[1], seqs[2]))
                                    for s in seqs:
                                        if s not in count_by_sequence:
                                            count_by_sequence[s] = 1
                                        else:
                                            count_by_sequence[s] += 1
                                        
                                triangle_nodes_all.add (triad)
                                
                            
                            
            
            
        self.clear_adjacency()
        print ("Found %d (%d) triangles in the network" % (len (triangles), len (triangle_nodes_all)), file = sys.stderr)
    
        return [ (t[0], t[1], t[2], sum([count_by_sequence [t[0]],count_by_sequence [t[1]],count_by_sequence [t[2]]]))  for t in triangles] 
                            

    def test_edge_support (self, sequence_file_name, triangles, hy_instance = None, p_value_cutoff = 0.05):

    
        evaluator = partial (_test_edge_support, sequence_file_name = sequence_file_name, hy_instance = hy_instance, p_value_cutoff = p_value_cutoff)
        #processed_objects = evaluator (triangles)
        
        blocked = [triangles[k : k + 256] for k in range (0, len (triangles), 256)]
        
        pool = multiprocessing.Pool()
        #print (multiprocessing.cpu_count())
        processed_objects = pool.map(evaluator, blocked)
        pool.close ()
        pool.join ()
        
        seqs_to_edge = {}
        for e in self.edges:
            if e.sequences:
                seqs_to_edge [e.sequences] = e
                
                
        processed_objects = sorted ([k for p in processed_objects for k in p], key = lambda x: x[0][3])   
                
        edges_removed = set()
        must_keep     = set()     
        
              

        for t, p_values in processed_objects:
            seq_id = t[:3]
            edges = [None,None,None]
            for pair_index,pair in enumerate(((0,1),(0,2),(1,2))):
                for seq_tag in [(seq_id[pair[0]], seq_id[pair[1]]), (seq_id[pair[1]], seq_id[pair[0]])]:
                    if seq_tag in seqs_to_edge:
                        edges[pair_index] = seqs_to_edge[seq_tag]
                        break
        
            for i in range (3):
               edges[i].edge_reject_p = max (p_values[2-i], edges[i].edge_reject_p)
           
            max_p = max (p_values)
            if max_p > p_value_cutoff:
                remove_index = 2 - p_values.index(max_p)
                if edges[remove_index] not in edges_removed and edges[remove_index] not in must_keep:
                    edges[remove_index].is_unsupported = True
                    edges_removed.add (edges[remove_index])
                    for e in [k for k in range (3) if k!=remove_index]:
                        must_keep.add (edges[e])
                            
            
        
                
    
    def fit_degree_distribution (self, degree_option = None, hy_instance = None):
        if hy_instance is None:
            hy_instance = hy.HyphyInterface ();
        script_path = os.path.realpath(__file__)
        hbl_path =  os.path.join(os.path.dirname(script_path), "data", "HBL", "DegreeDistributions.bf")
        if degree_option == 'indegree':
            all_deg = self.get_degree_distribution(indegree=True)
        elif degree_option == 'outdegree':
            all_deg = self.get_degree_distribution(outdegree=True)
        else:
            if degree_option is None:
                all_deg = self.get_degree_distribution()
            else:
                all_deg = degree_option

        hy_instance.queuevar ('allDegs', all_deg)
        hy_instance.runqueue (batchfile = hbl_path)
        bestDistro = hy_instance.getvar ('BestDistro',hy.HyphyInterface.STRING)
        rho = {}
        bic = {}
        p = {}
        fitted = {}
        rho_ci = {}
        for name in ('Waring', 'Yule', 'Pareto', 'Negative Binomial'):
            try:
                rho[name] = hy_instance.getvar (name,hy.HyphyInterface.NUMBER)
            except:
                rho[name] = None
            try:
                rho_ci[name] = hy_instance.getvar (name + "_rho_ci",hy.HyphyInterface.MATRIX)
            except:
                rho_ci[name] = None
            try:
                bic[name] = hy_instance.getvar (name + "_BIC",hy.HyphyInterface.NUMBER)
            except:
                bic[name] = None
            try:
                p[name] = hy_instance.getvar (name + "_p",hy.HyphyInterface.NUMBER)
            except:
                p[name] = None
            try:
                fitted[name] = hy_instance.getvar (name + "_PDF",hy.HyphyInterface.MATRIX)
            except:
                fitted[name] = None
        return {'Best' : bestDistro, 'rho': rho, 'BIC' : bic, 'p' : p, 'fitted': fitted, 'degrees': all_deg, 'rho_ci': rho_ci}

    def simulate_treatment (self,treated_nodes, node_neighs_out, node_neighs_in, node_neighs_undirected, removal_rate = 1.):
        removed_nodes = set()
        handled_nodes = set()
        new_nodes     = set(treated_nodes)

        while True:
            for a_removed_node in removed_nodes:
                if a_removed_node not in handled_nodes:
                    nbhd = (node_neighs_out[a_removed_node] - removed_nodes - new_nodes).union(node_neighs_undirected[a_removed_node])
                    if (len (nbhd)) :
                        sel = list (nbhd)
                        random.shuffle(sel)
                        for nb_node in sel:
                            potential_sources = node_neighs_in[nb_node] - removed_nodes - new_nodes
                            local_undir       = node_neighs_undirected[nb_node] - removed_nodes - new_nodes

                            potential_sources.add (a_removed_node)
                            pure_links = len (potential_sources)
                            potential_sources = potential_sources.union(local_undir)
                            src_count = 0.5*(len(potential_sources)+pure_links)

                            if random.random() <= (0.5 if nb_node in node_neighs_undirected[a_removed_node] else 1.)/(src_count):
                                if nb_node not in treated_nodes or random.random () < removal_rate:
                                    new_nodes.add (nb_node)

                    handled_nodes.add (a_removed_node)

            if len (new_nodes) == 0:
                break
            else:
                removed_nodes.update(new_nodes)
            new_nodes = set ()

        removed_nodes = removed_nodes.difference (treated_nodes)
        return removed_nodes

    def get_degree_distribution (self, **kwargs):
        degree_distribution = []

        subset = None
        if 'subset' in kwargs:
            subset = kwargs ['subset']

        directed = False
        if 'directed' in kwargs:
            directed = bool (kwargs ['directed'])

        outdegree = False
        if 'outdegree' in kwargs:
            outdegree = bool (kwargs ['outdegree'])

        indegree = False
        if 'indegree' in kwargs:
            indegree = bool (kwargs ['indegree'])


        per_year_fu = None
        if 'peryear' in kwargs:
            per_year_fu = int (kwargs ['peryear'])

        per_node    = None
        if 'storenodes' in kwargs:
            per_node = kwargs['storenodes']

        if self.adjacency_list == None or ((directed or outdegree or indegree) and self.type_of_adjacency_list () == 'patient') or (not (directed or outdegree or indegree) and self.type_of_adjacency_list () == 'edge'):
            #print 'Redo'
            self.compute_adjacency(directed or outdegree or indegree)

        max_diff    = None
        if 'max_diff' in kwargs and 'directed':
            if isinstance (kwargs['max_diff'], int):
                max_diff = datetime.timedelta (days = int (kwargs['max_diff']))

        for node in self.adjacency_list:
            if subset and node not in subset:
                continue

            if directed or outdegree or indegree:
                this_degree = 0
                for an_edge in self.adjacency_list[node]:
                    dir = an_edge.compute_direction()

                    connect_me = False
                    if outdegree:
                        connect_me = dir is not None and dir == node
                    elif indegree:
                        connect_me = dir is not None and dir != node
                    else:
                        connect_me = dir is not None and dir == node

                    if connect_me:
                        if max_diff:
                            diff = an_edge.chrono_length_days()
                            if diff == None or diff <= max_diff:
                                this_degree += 1
                        else:
                            this_degree += 1
            else:
                this_degree = len (self.adjacency_list[node])

            if per_year_fu:
                degree_distribution.append (this_degree/float(per_year_fu - node.get_baseline_date () + 1))
            else:
                if len (degree_distribution) < this_degree + 1:
                    for k in range (len (degree_distribution), this_degree+1):
                        degree_distribution.append (0)
                degree_distribution[this_degree] += 1

            if per_node is not None:
                per_node [node] = this_degree

        #print "Degree : %s " % str (degree_distribution)

        if 'transform' in kwargs and not per_year_fu:
            if kwargs['transform'] == 'NetworkStat':
                deg = 0
                for i,d in enumerate(degree_distribution):
                    deg += (i+1)*d

                return deg

            normalizer = 1./sum (degree_distribution)
            degree_distribution = [k * normalizer for k in degree_distribution]

            if kwargs['transform'] == 'CDF' or kwargs['transform'] == 'LogCDF':
                cdf = copy(degree_distribution)
                for k in range (1, len (degree_distribution)):
                    cdf[k] += cdf[k-1]
                degree_distribution = [cdf[k] - degree_distribution[k] for k in range (len (degree_distribution))]
                if kwargs['transform'] == 'LogCDF':
                    return [log (1-k) for k in degree_distribution [1:]]

        if subset:
            if per_year_fu:
                for k in range (len (subset) - len (degree_distribution)):
                    degree_distribution.append (0.0)
            else:
                degree_distribution[0] += len (subset) - sum (degree_distribution)
            return degree_distribution

        return degree_distribution [1:]
