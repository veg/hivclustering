#!/usr/bin/env python3.2

from operator import itemgetter
import json
import csv
import os.path
import nose


from hivclustering import *
#from networkbuild import *

def test_centrality():
    ''' Creates a Krackhardt kite and ensures betweenness is correct '''

    #Create a transmission network
    network = transmission_network()

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

    centralities = {}
    [centralities.update(network.betweenness_centrality(i)) for i in range(10)]
    assert len(centralities.keys()) == 10
    assert centralities['Heather'] - .388888 < .0001

