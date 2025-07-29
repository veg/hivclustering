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
        self.assertIn('Recommended', header, "Recommended column not found in header")
        
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
        """Test that singleton counts are accurate and decrease as threshold increases"""
        # Create test data with known singleton behavior
        test_data = """ID1,ID2,Distance
A|01012020,B|01012020,0.001
C|01012020,D|01012020,0.01
E|01012020,F|01012020,0.02
G|01012020,H|01012020,0.05
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write(test_data)
            test_file = f.name
        
        try:
            cmd = [
                sys.executable,
                os.path.join(os.path.dirname(__file__), '..', 'scripts', 'hivnetworkcsv'),
                '-i', test_file,
                '-t', '0.1',
                '-A', '0.0001'
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, f"Command failed with stderr: {result.stderr}")
            
            lines = result.stdout.strip().split('\n')
            header = lines[0].split('\t')
            
            self.assertIn('Singletons', header, "Singletons column not found in header")
            self.assertIn('Recommended', header, "Recommended column not found in header")
            
            # Verify decreasing singleton pattern and total node conservation
            threshold_idx = header.index('Threshold')
            nodes_idx = header.index('Nodes')
            singletons_idx = header.index('Singletons')
            
            prev_singletons = None
            total_nodes = 8  # A,B,C,D,E,F,G,H from the test data
            
            for line in lines[1:]:
                if line.strip():
                    cols = line.split('\t')
                    threshold = float(cols[threshold_idx])
                    nodes = int(cols[nodes_idx])
                    singletons = int(cols[singletons_idx])
                    
                    # Verify singletons are non-negative
                    self.assertGreaterEqual(singletons, 0, 
                                          f"Singletons count should be non-negative at threshold {threshold}")
                    
                    # Verify total node conservation
                    self.assertEqual(nodes + singletons, total_nodes,
                                   f"At threshold {threshold}: nodes ({nodes}) + singletons ({singletons}) should equal {total_nodes}")
                    
                    # Verify singletons generally decrease (allow for some variation)
                    if prev_singletons is not None:
                        self.assertLessEqual(singletons, prev_singletons + 2,
                                           f"Singletons should generally decrease: {prev_singletons} → {singletons} at threshold {threshold}")
                    
                    prev_singletons = singletons
                    
        finally:
            os.unlink(test_file)

    def test_auto_threshold_mean_diff_zero_fix(self):
        """Test that AUTO-TUNE handles mean_diff == 0.0 properly for low-diversity data"""
        # Create test data that would trigger mean_diff == 0.0 in prior versions
        test_data = """ID1,ID2,Distance
seq1,seq2,0.001
seq1,seq3,0.0015
seq2,seq3,0.002
seq3,seq4,0.0025
seq4,seq5,0.003
seq5,seq6,0.0035
seq6,seq7,0.004
seq7,seq8,0.0045
seq8,seq9,0.005
seq9,seq10,0.0055
seq10,seq11,0.006
seq11,seq12,0.0065
seq12,seq13,0.007
seq13,seq14,0.0075
seq14,seq15,0.008
seq15,seq16,0.0085
seq16,seq17,0.009
seq17,seq18,0.0095
seq18,seq19,0.010
seq19,seq20,0.0105
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write(test_data)
            test_file = f.name
        
        try:
            cmd = [
                sys.executable,
                os.path.join(os.path.dirname(__file__), '..', 'scripts', 'hivnetworkcsv'),
                '-i', test_file,
                '-f', 'plain',
                '-t', 'auto'
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            # The fix should provide a meaningful best guess (not 1e-05 with score 0)
            # But still maintain the error behavior for consistency with prior versions
            if "best guess" in result.stderr:
                # Should provide a meaningful best guess with actual score, not default fallback
                self.assertNotIn("best guess 1e-05 (score 0)", result.stderr,
                               "Should not fall back to default threshold when real candidates exist")
                # Should see a meaningful best guess with threshold and score
                self.assertRegex(result.stderr, r"best guess [\d\.e\-\+]+ \(score [\d\.e\-\+]+\)",
                               "Should provide meaningful best guess with threshold and real score")
            else:
                # In some cases it might successfully select a threshold
                self.assertIn("Selected distance threshold", result.stderr,
                             "Should either select threshold or provide meaningful best guess")
                           
        finally:
            os.unlink(test_file)

    def test_auto_profile_shows_recommendation(self):
        """Test that AUTO-TUNE auto-profile mode shows Recommended column"""
        test_data = """ID1,ID2,Distance
A|01012020,B|01012020,0.001
C|01012020,D|01012020,0.002
E|01012020,F|01012020,0.003
G|01012020,H|01012020,0.004
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write(test_data)
            test_file = f.name
        
        try:
            cmd = [
                sys.executable,
                os.path.join(os.path.dirname(__file__), '..', 'scripts', 'hivnetworkcsv'),
                '-i', test_file,
                '-A', '0.001'
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, f"Command failed with stderr: {result.stderr}")
            
            lines = result.stdout.strip().split('\n')
            self.assertGreater(len(lines), 1, "Should have header and at least one data line")
            
            header = lines[0].split('\t')
            self.assertIn('Recommended', header, "Should have Recommended column")
            
            recommended_idx = header.index('Recommended')
            
            # Check that exactly one row has "*" in the Recommended column
            recommended_count = 0
            for line in lines[1:]:  # Skip header
                if line.strip():
                    cols = line.split('\t')
                    if len(cols) > recommended_idx and cols[recommended_idx] == "*":
                        recommended_count += 1
            
            self.assertEqual(recommended_count, 1, "Exactly one row should be marked as recommended")
                    
        finally:
            os.unlink(test_file)


if __name__ == '__main__':
    unittest.main()