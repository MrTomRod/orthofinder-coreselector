"""
Simple test for OrthoFinder CoreSelector functionality.
"""

import os
import unittest

from orthofinder_coreselector import select_core_set


class TestOrthoFinderCoreSelector(unittest.TestCase):
    """Test OrthoFinder CoreSelector functionality."""

    def test_reduced_selection(self):
        """Test basic genome selection functionality."""

        os.makedirs("tests/data/out/reduced", exist_ok=True)
        target_core_size = 20

        # Run the core selection
        selected_genomes = select_core_set(
            metadata_file="tests/data/genomes.tsv",
            distance_matrix="tests/data/ani-distance-matrix.csv",
            target_core_size=target_core_size,
            output_file="tests/data/out/reduced/selected_genomes.txt",
            visualization_file="tests/data/out/reduced/visualisation.svg",
            config_file="tests/data/config.json",
            max_genomes=200
        )

        # Basic assertions
        self.assertIsInstance(selected_genomes, list)
        self.assertEqual(len(selected_genomes), target_core_size)
        self.assertTrue(all(isinstance(genome, str) for genome in selected_genomes))

    def test_reduced_selection_no_metadata(self):
        """Test basic genome selection functionality."""

        os.makedirs("tests/data/out/reduced", exist_ok=True)
        target_core_size = 20

        # Run the core selection
        selected_genomes = select_core_set(
            distance_matrix="tests/data/ani-distance-matrix.csv",
            target_core_size=target_core_size,
            output_file="tests/data/out/reduced/selected_genomes_no_metadata.txt",
            visualization_file="tests/data/out/reduced/visualisation_no_metadata.svg",
            max_genomes=200
        )

        # Basic assertions
        self.assertIsInstance(selected_genomes, list)
        self.assertEqual(len(selected_genomes), target_core_size)
        self.assertTrue(all(isinstance(genome, str) for genome in selected_genomes))

if __name__ == "__main__":
    unittest.main()
