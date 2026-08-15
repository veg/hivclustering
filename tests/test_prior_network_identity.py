#!/usr/bin/env python3

"""
    Tests for how hivnetworkcsv matches the current network against a prior network (-P),
    i.e. which nodes are reported as new, removed, or moved between clusters.

    Sequences from the same person are named PERSON|SEQUENCE..., so replacing a person's
    sequence must not make that person a new node (--prior-identity entity, the default),
    while --prior-identity sequence preserves the pre-1.9.10 behavior of matching on the
    full sequence ID.
"""

import unittest
import sys
import os
import json
import shutil
import subprocess
import tempfile

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from hivclustering import entity_id_from_string

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
HIVNETWORKCSV = os.path.join(REPO_ROOT, 'scripts', 'hivnetworkcsv')

# a three person cluster (A, B, C), a two person cluster (D, E)
PRIOR_LINKS = [
    ("A|s1", "B|s1", 0.001),
    ("B|s1", "C|s1", 0.002),
    ("D|s1", "E|s1", 0.001),
]


def write_links(path, links):
    with open(path, 'w') as fh:
        fh.write("ID1,ID2,Distance\n")
        for id1, id2, distance in links:
            fh.write("%s,%s,%g\n" % (id1, id2, distance))
    return path


class TestPriorNetworkIdentity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.workdir = tempfile.mkdtemp()
        cls.prior_json = cls.build_network(write_links(os.path.join(cls.workdir, 'prior.csv'), PRIOR_LINKS),
                                           os.path.join(cls.workdir, 'prior.json'))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.workdir, ignore_errors=True)

    @classmethod
    def build_network(cls, links_csv, output_json, prior_json=None, identity=None):
        cmd = [sys.executable, HIVNETWORKCSV, '-i', links_csv, '-f', 'plain', '-j',
               '-t', '0.015', '--no-degree-fit', '-O', output_json]
        if prior_json:
            cmd += ['-P', prior_json]
        if identity:
            cmd += ['--prior-identity', identity]

        env = dict(os.environ)
        # make sure the working copy, not an installed release, is what gets tested
        env['PYTHONPATH'] = REPO_ROOT + os.pathsep + env.get('PYTHONPATH', '')

        result = subprocess.run(cmd, capture_output=True, text=True, env=env)
        assert result.returncode == 0, "hivnetworkcsv failed with stderr: %s" % result.stderr
        return output_json

    def compare_to_prior(self, name, links, identity=None):
        """Build a network from links, compare it to the shared prior network, return its JSON"""
        links_csv = write_links(os.path.join(self.workdir, '%s.csv' % name), links)
        output_json = os.path.join(self.workdir, '%s-%s.json' % (name, identity or 'default'))
        self.build_network(links_csv, output_json, prior_json=self.prior_json, identity=identity)
        with open(output_json) as fh:
            return json.load(fh)

    @staticmethod
    def new_nodes(network):
        return set([n['id'] for n in network['Nodes'] if 'new_node' in n.get('attributes', [])])

    @staticmethod
    def cluster_of(network, node_id):
        return [n['cluster'] for n in network['Nodes'] if n['id'] == node_id][0]

    def test_entity_id_from_string(self):
        self.assertEqual(entity_id_from_string("GA00S000163429-5|GA00I000181113-7~3~PR_RT~1212"), "GA00S000163429-5")
        self.assertEqual(entity_id_from_string("AL00S000501064-0"), "AL00S000501064-0")
        self.assertEqual(entity_id_from_string("person|seq|extra"), "person")

    def test_resequenced_person_is_not_new(self):
        """A person whose sequence ID changed is the same person, not a new node"""
        resequenced = [("A|s2", "B|s1", 0.001)] + PRIOR_LINKS[1:]
        network = self.compare_to_prior('resequenced', resequenced)

        self.assertEqual(self.new_nodes(network), set(),
                         "re-sequencing an existing person must not flag them as a new node")
        self.assertNotIn('Notes', network,
                         "re-sequencing an existing person must not report a removed node")
        self.assertEqual(self.cluster_of(network, 'A|s2'), self.cluster_of(network, 'B|s1'))

        cluster_description = network['Cluster description'][str(self.cluster_of(network, 'A|s2'))]
        self.assertNotIn('new_nodes', cluster_description,
                         "re-sequencing an existing person must not count as cluster growth")

        edge = [e for e in network['Edges'] if set(e['sequences']) == set(['A|s2', 'B|s1'])][0]
        self.assertNotIn('added-to-prior', edge['attributes'],
                         "a link that existed before must not be reported as added")

    def test_resequenced_person_is_new_under_sequence_identity(self):
        """--prior-identity sequence keeps the pre-1.9.10 (full sequence ID) behavior"""
        resequenced = [("A|s2", "B|s1", 0.001)] + PRIOR_LINKS[1:]
        network = self.compare_to_prior('resequenced', resequenced, identity='sequence')

        self.assertEqual(self.new_nodes(network), set(['A|s2']))
        self.assertIn('Notes', network)
        self.assertTrue(any('Removed 1 nodes' in note for note in network['Notes']))

    def test_additional_sequence_for_existing_person_is_not_new(self):
        """An existing person gaining another specimen is not network growth"""
        network = self.compare_to_prior('additional', PRIOR_LINKS + [("A|s1", "A|s2", 0.0001)])

        self.assertEqual(self.new_nodes(network), set())
        self.assertEqual(self.cluster_of(network, 'A|s2'), self.cluster_of(network, 'A|s1'),
                         "an added sequence should join the cluster of the person it belongs to")

    def test_genuinely_new_person_is_new(self):
        """A person who was not in the prior network is still flagged, under either identity"""
        with_new_person = PRIOR_LINKS + [("D|s1", "F|s1", 0.001)]

        for identity in (None, 'sequence'):
            network = self.compare_to_prior('newperson', with_new_person, identity=identity)
            self.assertEqual(self.new_nodes(network), set(['F|s1']))

            cluster_description = network['Cluster description'][str(self.cluster_of(network, 'F|s1'))]
            self.assertEqual(cluster_description['new_nodes'], 1)

            edge = [e for e in network['Edges'] if set(e['sequences']) == set(['D|s1', 'F|s1'])][0]
            self.assertIn('added-to-prior', edge['attributes'])

    def test_departed_person_is_still_reported_as_removed(self):
        """Dropping a person entirely is real attrition and must still be reported"""
        network = self.compare_to_prior('departed', [l for l in PRIOR_LINKS if 'D|s1' not in l])

        self.assertIn('Notes', network)
        self.assertTrue(any('Removed 2 nodes' in note for note in network['Notes']),
                        "D and E both left the network")

    def test_ids_without_an_entity_prefix_are_unaffected(self):
        """Networks whose node IDs carry no entity prefix behave identically either way"""
        links = [(id1.replace('|', '_'), id2.replace('|', '_'), d) for id1, id2, d in PRIOR_LINKS]
        prior = self.build_network(write_links(os.path.join(self.workdir, 'flat-prior.csv'), links),
                                   os.path.join(self.workdir, 'flat-prior.json'))

        grown = links + [("D_s1", "F_s1", 0.001)]
        links_csv = write_links(os.path.join(self.workdir, 'flat.csv'), grown)

        networks = {}
        for identity in ('entity', 'sequence'):
            output_json = os.path.join(self.workdir, 'flat-%s.json' % identity)
            self.build_network(links_csv, output_json, prior_json=prior, identity=identity)
            with open(output_json) as fh:
                networks[identity] = json.load(fh)

        self.assertEqual(self.new_nodes(networks['entity']), set(['F_s1']))
        self.assertEqual(networks['entity']['Nodes'], networks['sequence']['Nodes'])
        self.assertEqual(networks['entity']['Cluster description'], networks['sequence']['Cluster description'])


if __name__ == '__main__':
    unittest.main()
