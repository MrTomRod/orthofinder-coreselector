"""
Tests for data loading functionality.
"""

# Add src to path for imports
import sys
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))

from orthofinder_coreselector.data_loading import load_distance_matrix, load_metadata


class TestDataLoading(unittest.TestCase):
    """Test data loading functionality."""

    def test_load_metadata(self):
        """Test metadata loading."""
        metadata_file = Path(__file__).parent / "data" / "genomes.tsv"
        assert  metadata_file.exists()

        metadata = load_metadata(str(metadata_file), index_col="Identifier")

        # Basic assertions
        self.assertIsInstance(metadata, pd.DataFrame)
        self.assertGreater(len(metadata), 0)
        # Identifier becomes the index, so check index name instead
        self.assertEqual(metadata.index.name, 'Identifier')
        # Check that we have some expected columns (at least one should exist)
        expected_cols = ['Representative', 'Assembly Size', 'Assembly Nr Scaffolds', 'BUSCO [S]', 'Species']
        found_cols = [col for col in expected_cols if col in metadata.columns]
        self.assertGreater(len(found_cols), 0, \
            f"Expected at least one of {expected_cols} but found: {list(metadata.columns)}")

    def test_load_distance_matrix(self):
        """Test distance matrix loading."""
        distance_file = Path(__file__).parent / "data" / "ani-distance-matrix.csv"

        if distance_file.exists():
            distance_data = load_distance_matrix(str(distance_file))

            # Basic assertions
            self.assertIsInstance(distance_data, pd.DataFrame)
            self.assertGreater(len(distance_data), 0)
            self.assertEqual(len(distance_data.index), len(distance_data.columns))
            # Check that diagonal is zero (self-distance)
            for i in range(len(distance_data)):
                self.assertEqual(distance_data.iloc[i, i], 0.0)


if __name__ == "__main__":
    unittest.main()
