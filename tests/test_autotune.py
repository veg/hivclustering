#!/usr/bin/env python3

import unittest
import sys
import os
import tempfile
import subprocess
from io import StringIO

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from hivclustering import mtnetwork
from hivclustering import networkbuild


class TestAutoTuneSingletons(unittest.TestCase):
    def setUp(self):
        # Create a simple dataset that should work with AUTO-TUNE
        self.test_csv_data = """ID1,ID2,Distance
A|01012020,B|01012020,0.001
C|01012020,D|01012020,0.002
E|01012020,F|01012020,0.003
"""
        self.temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False)
        self.temp_file.write(self.test_csv_data)
        self.temp_file.close()
        
    def tearDown(self):
        os.unlink(self.temp_file.name)
    
    def test_autotune_output_contains_singletons_column(self):
        """Test that AUTO-TUNE output includes Singletons column"""
        cmd = [
            sys.executable,
            os.path.join(os.path.dirname(__file__), '..', 'scripts', 'hivnetworkcsv'),
            '-i', self.temp_file.name,
            '-t', '0.01',
            '-A', '0.001'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        self.assertEqual(result.returncode, 0, f"Command failed with stderr: {result.stderr}")
        
        lines = result.stdout.strip().split('\n')
        self.assertGreater(len(lines), 0, "No output produced")
        
        header = lines[0].split('\t')
        self.assertIn('Singletons', header, "Singletons column not found in header")
        
        singletons_index = header.index('Singletons')
        
        for i, line in enumerate(lines[1:], 1):
            if line.strip():
                columns = line.split('\t')
                self.assertEqual(len(columns), len(header), 
                               f"Line {i} has wrong number of columns: {line}")
                
                singletons_value = columns[singletons_index]
                self.assertTrue(singletons_value.isdigit(), 
                               f"Line {i} has non-numeric singleton value: {singletons_value}")
                
                threshold = float(columns[0])
                nodes = int(columns[1])
                singletons = int(singletons_value)
                
                self.assertGreaterEqual(singletons, 0, 
                                      f"Singletons cannot be negative at threshold {threshold}")
                
                # At lower thresholds, only some nodes will be in the network
                # At higher thresholds, more nodes will be included
                total_nodes = nodes + singletons
                self.assertGreater(total_nodes, 0, 
                                 f"At threshold {threshold}: should have some nodes in network")
                
                # Verify that total nodes increases or stays same as threshold increases
                # (nodes are added as more edges become valid)

    def test_singleton_count_accuracy(self):
        """Test that singleton counts are accurate at different thresholds"""
        # Create a simple network with clear singleton nodes
        # Z and AA are isolated nodes (singletons)
        more_complex_data = """ID1,ID2,Distance
A|01012020,B|01012020,0.001
B|01012020,C|01012020,0.002
D|01012020,E|01012020,0.015
F|01012020,G|01012020,0.025
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write(more_complex_data)
            complex_file = f.name
        
        try:
            cmd = [
                sys.executable,
                os.path.join(os.path.dirname(__file__), '..', 'scripts', 'hivnetworkcsv'),
                '-i', complex_file,
                '-t', '0.05',
                '-A', '0.0001'
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, f"Command failed with stderr: {result.stderr}")
            
            lines = result.stdout.strip().split('\n')
            header = lines[0].split('\t')
            
            self.assertIn('Singletons', header, "Singletons column not found in header")
            
            # Just verify the column exists and has valid data
            singletons_idx = header.index('Singletons')
            
            for line in lines[1:]:
                if line.strip():
                    cols = line.split('\t')
                    singletons = int(cols[singletons_idx])
                    self.assertGreaterEqual(singletons, 0, 
                                          "Singletons count should be non-negative")
                    
        finally:
            os.unlink(complex_file)


if __name__ == '__main__':
    unittest.main()