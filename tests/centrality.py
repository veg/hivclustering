#!/usr/bin/env python3

from operator import itemgetter
import json
import csv
import os.path
import nose
import functools


from hivclustering import *
#from networkbuild import *
network = transmission_network()

def setup():
    ''' Creates a Krackhardt kite and ensures betweenness is correct '''
    #Create a transmission network

    #Add the nodes
    attrib = None
    network.insert_patient('Andre',False, False, attrib)
    network.insert_patient('Beverly',False, False, attrib)
    network.insert_patient('Carol',False, False, attrib)
    network.insert_patient('Diane',False, False, attrib)
    network.insert_patient('Ed',False, False, attrib)
    network.insert_patient('Fernando',False, False, attrib)
    network.insert_patient('Garth',False, False, attrib)
    network.insert_patient('Heather',False, False, attrib)
    network.insert_patient('Ike',False, False, attrib)

    #Add the edges
    network.add_an_edge('Andre', 'Beverly', 1, parsePlain)
    network.add_an_edge('Andre', 'Carol', 1, parsePlain)
    network.add_an_edge('Andre', 'Diane', 1, parsePlain)
    network.add_an_edge('Andre', 'Fernando', 1, parsePlain)
    network.add_an_edge('Beverly', 'Diane', 1, parsePlain)
    network.add_an_edge('Beverly', 'Ed', 1, parsePlain)
    network.add_an_edge('Beverly', 'Garth', 1, parsePlain)
    network.add_an_edge('Carol', 'Diane', 1, parsePlain)
    network.add_an_edge('Carol', 'Fernando', 1, parsePlain)
    network.add_an_edge('Diane', 'Ed', 1, parsePlain)
    network.add_an_edge('Diane', 'Fernando', 1, parsePlain)
    network.add_an_edge('Diane', 'Garth', 1, parsePlain)
    network.add_an_edge('Ed', 'Garth', 1, parsePlain)
    network.add_an_edge('Fernando', 'Garth', 1, parsePlain)
    network.add_an_edge('Fernando', 'Heather', 1, parsePlain)
    network.add_an_edge('Garth', 'Heather', 1, parsePlain)
    network.add_an_edge('Heather', 'Ike', 1, parsePlain)
    network.add_an_edge('Ike', 'Jane', 1, parsePlain)


@nose.with_setup(setup=setup)
def test_centrality():
    ''' Ensure betweenness is correct '''
    assert network.betweenness_centrality('Heather') - .388888 < .0001

@nose.with_setup(setup=setup)
def test_degrees():
    ''' Ensure degrees are correct '''
    patients = [k for k,v in network.nodes.items()]
    patients = [(patient.id, patient.degree) for patient in patients]
    expected = [('Carol', 9), ('Ed', 9), ('Diane', 18), ('Jane', 3), ('Fernando', 15), ('Andre', 12), ('Ike', 6), ('Beverly', 12), ('Heather', 9), ('Garth', 15)]
    assert set(patients) == set(expected)

