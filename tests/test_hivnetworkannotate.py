#!/usr/bin/env python3

import unittest
import sys
import os
import tempfile
import subprocess
import json

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from hivclustering import networkbuild


class TestHivnetworkannotateTraceResults(unittest.TestCase):
    def setUp(self):
        # Create minimal network JSON with trace_results wrapper and version info
        self.network_with_trace_results = {
            "hivtrace_version": "0.10.1",
            "trace_results": {
                "Nodes": [
                    {"id": "patient1~seq1"},
                    {"id": "patient2~seq2"},
                    {"id": "patient3~seq3"}
                ],
                "Edges": [
                    {
                        "sequences": ["patient1~seq1", "patient2~seq2"],
                        "length": 0.005
                    }
                ],
                "Settings": {
                    "threshold": 0.015
                }
            },
            "other_metadata": "preserved_value"
        }

        # Create minimal network JSON without trace_results wrapper
        self.network_without_trace_results = {
            "Nodes": [
                {"id": "patient1~seq1"},
                {"id": "patient2~seq2"},
                {"id": "patient3~seq3"}
            ],
            "Edges": [
                {
                    "sequences": ["patient1~seq1", "patient2~seq2"],
                    "length": 0.005
                }
            ],
            "Settings": {
                "threshold": 0.015
            }
        }

        # Create minimal attributes data
        self.attributes_data = [
            {
                "patient1": "patient1",
                "test_field": "value1"
            },
            {
                "patient2": "patient2",
                "test_field": "value2"
            },
            {
                "patient3": "patient3",
                "test_field": "value3"
            }
        ]

        # Create fields file for annotation
        self.fields_data = {
            "test_field": {
                "label": "Test Field",
                "type": "String"
            },
            "keying": {
                "fields": ["patient1"],
                "delimiter": "~"
            }
        }

    def test_trace_results_wrapper_preserves_version(self):
        """Test that trace_results wrapper preserves hivtrace_version and other top-level keys"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as network_file:
            json.dump(self.network_with_trace_results, network_file)
            network_path = network_file.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as attr_file:
            json.dump(self.attributes_data, attr_file)
            attr_path = attr_file.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as fields_file:
            json.dump(self.fields_data, fields_file)
            fields_path = fields_file.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as output_file:
            output_path = output_file.name

        try:
            cmd = [
                sys.executable,
                os.path.join(os.path.dirname(__file__), '..', 'scripts', 'hivnetworkannotate'),
                '-n', network_path,
                '-a', attr_path,
                '-g', fields_path,
                '-o', output_path
            ]

            result = subprocess.run(cmd, capture_output=True, text=True)

            self.assertEqual(result.returncode, 0,
                           f"Command failed with stderr: {result.stderr}")

            # Read the output and verify structure
            with open(output_path, 'r') as f:
                output_json = json.load(f)

            # Verify top-level keys are preserved
            self.assertIn('hivtrace_version', output_json,
                         "hivtrace_version should be preserved in output")
            self.assertEqual(output_json['hivtrace_version'], "0.10.1",
                           "hivtrace_version value should match input")

            self.assertIn('other_metadata', output_json,
                         "other_metadata should be preserved in output")
            self.assertEqual(output_json['other_metadata'], "preserved_value",
                           "other_metadata value should match input")

            # Verify trace_results is still present
            self.assertIn('trace_results', output_json,
                         "trace_results should be present in output")

            # Verify trace_results contains the annotated network
            self.assertIn('Nodes', output_json['trace_results'],
                         "trace_results should contain Nodes")
            self.assertIn('Edges', output_json['trace_results'],
                         "trace_results should contain Edges")

            # Verify patient_attribute_schema was added
            self.assertIn('patient_attribute_schema', output_json['trace_results'],
                         "patient_attribute_schema should be added to trace_results")

        finally:
            os.unlink(network_path)
            os.unlink(attr_path)
            os.unlink(fields_path)
            os.unlink(output_path)

    def test_network_without_trace_results_unchanged(self):
        """Test that networks without trace_results wrapper work correctly"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as network_file:
            json.dump(self.network_without_trace_results, network_file)
            network_path = network_file.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as attr_file:
            json.dump(self.attributes_data, attr_file)
            attr_path = attr_file.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as fields_file:
            json.dump(self.fields_data, fields_file)
            fields_path = fields_file.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as output_file:
            output_path = output_file.name

        try:
            cmd = [
                sys.executable,
                os.path.join(os.path.dirname(__file__), '..', 'scripts', 'hivnetworkannotate'),
                '-n', network_path,
                '-a', attr_path,
                '-g', fields_path,
                '-o', output_path
            ]

            result = subprocess.run(cmd, capture_output=True, text=True)

            self.assertEqual(result.returncode, 0,
                           f"Command failed with stderr: {result.stderr}")

            # Read the output and verify structure
            with open(output_path, 'r') as f:
                output_json = json.load(f)

            # Verify NO trace_results wrapper (flat structure)
            self.assertNotIn('trace_results', output_json,
                           "Output should not have trace_results wrapper")

            # Verify network data is at top level
            self.assertIn('Nodes', output_json,
                         "Nodes should be at top level")
            self.assertIn('Edges', output_json,
                         "Edges should be at top level")

            # Verify patient_attribute_schema was added
            self.assertIn('patient_attribute_schema', output_json,
                         "patient_attribute_schema should be added at top level")

        finally:
            os.unlink(network_path)
            os.unlink(attr_path)
            os.unlink(fields_path)
            os.unlink(output_path)


if __name__ == '__main__':
    unittest.main()
