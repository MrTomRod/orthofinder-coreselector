"""
Tests for edge cases in plotting functionality.
"""

import logging

# Add src to path for imports
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))

from orthofinder_coreselector.plots import create_visualization


class TestPlotsEdgeCases(unittest.TestCase):
    """Test edge cases in plotting functionality."""

    def setUp(self):
        """Set up test data."""
        # Disable logging for tests
        logging.getLogger("orthofinder_coreselector").setLevel(logging.CRITICAL)

        # Create test metadata
        self.metadata = pd.DataFrame({
            'BUSCO [%S]': [95.0, 90.0, 85.0],
            'Assembly Nr Scaffolds': [100, 200, 300],
            'Assembly Size': [5000000, 6000000, 7000000],
            'Representative': [True, False, False],
            'Species': ['E. coli', 'S. aureus', 'B. subtilis']
        }, index=pd.Index(['genome1', 'genome2', 'genome3'], name='Identifier'))

        # Create test distance matrix
        self.distance_data = pd.DataFrame({
            'genome1': [0.0, 0.1, 0.2],
            'genome2': [0.1, 0.0, 0.15],
            'genome3': [0.2, 0.15, 0.0]
        }, index=['genome1', 'genome2', 'genome3'])

        # Test quality scores
        self.quality_scores = {
            'genome1': 1.5,
            'genome2': 0.8,
            'genome3': 0.3
        }

        # Test selected genomes
        self.selected_genomes = ['genome1', 'genome2']

        # Test config with multiple column types
        self.config = {
            'columns': {
                'BUSCO [%S]': {
                    'type': 'numeric',
                    'weight': 1.0,
                    'plot': True,
                    'display_unit': '%',
                    'decimal_places': 1
                },
                'Assembly Size': {
                    'type': 'numeric',
                    'weight': 0.5,
                    'plot': True,
                    'display_unit': 'bp',
                    'scale_factor': 1000000,  # Convert to Mbp
                    'decimal_places': 1
                },
                'Representative': {
                    'type': 'boolean',
                    'true_value': True,
                    'weight': 1.0,
                    'plot': True
                },
                'Species': {
                    'type': 'text',
                    'plot': True
                }
            }
        }

    def test_create_visualization_no_common_genomes(self):
        """Test visualization creation with no common genomes between metadata and distance matrix."""
        # Create distance matrix with different genome IDs
        bad_distance = pd.DataFrame({
            'other1': [0.0, 0.1, 0.2],
            'other2': [0.1, 0.0, 0.15],
            'other3': [0.2, 0.15, 0.0]
        }, index=['other1', 'other2', 'other3'])

        with tempfile.NamedTemporaryFile(suffix='.svg') as f:
            # Should raise ValueError with clear error message
            with self.assertRaises(ValueError) as cm:
                create_visualization(
                    metadata=self.metadata,
                    distance_data=bad_distance,
                    selected_genomes=self.selected_genomes,
                    quality_scores=self.quality_scores,
                    output_file=f.name,
                    config=self.config
                )
            self.assertIn("No common genomes found", str(cm.exception))

    def test_create_visualization_with_cluster_assignments(self):
        """Test visualization with cluster assignments."""
        cluster_assignments = {
            'genome1': 0,
            'genome2': 1,
            'genome3': 0
        }

        with tempfile.NamedTemporaryFile(suffix='.svg') as f:
            create_visualization(
                metadata=self.metadata,
                distance_data=self.distance_data,
                selected_genomes=self.selected_genomes,
                quality_scores=self.quality_scores,
                output_file=f.name,
                config=self.config,
                cluster_assignments=cluster_assignments
            )
            self.assertGreater(Path(f.name).stat().st_size, 1000)

    def test_create_visualization_with_filtered_genomes(self):
        """Test visualization with filtered genomes list."""
        filtered_genomes = ['genome1', 'genome3']  # genome2 was filtered out

        with tempfile.NamedTemporaryFile(suffix='.svg') as f:
            create_visualization(
                metadata=self.metadata,
                distance_data=self.distance_data,
                selected_genomes=['genome1'],  # Only genome1 selected from filtered
                quality_scores={'genome1': 1.5, 'genome3': 0.3},
                output_file=f.name,
                config=self.config,
                filtered_genomes=filtered_genomes
            )
            self.assertGreater(Path(f.name).stat().st_size, 1000)

    def test_create_visualization_missing_columns_in_metadata(self):
        """Test visualization when config references columns not in metadata."""
        config_with_missing = self.config.copy()
        config_with_missing['columns']['NonExistentColumn'] = {
            'type': 'numeric',
            'weight': 1.0,
            'plot': True
        }

        with tempfile.NamedTemporaryFile(suffix='.svg') as f:
            create_visualization(
                metadata=self.metadata,
                distance_data=self.distance_data,
                selected_genomes=self.selected_genomes,
                quality_scores=self.quality_scores,
                output_file=f.name,
                config=config_with_missing
            )
            self.assertGreater(Path(f.name).stat().st_size, 1000)

    def test_create_visualization_with_nan_values(self):
        """Test visualization with NaN values in metadata."""
        metadata_with_nan = self.metadata.copy()
        metadata_with_nan.loc['genome2', 'BUSCO [%S]'] = np.nan
        metadata_with_nan.loc['genome3', 'Assembly Size'] = np.nan

        with tempfile.NamedTemporaryFile(suffix='.svg') as f:
            create_visualization(
                metadata=metadata_with_nan,
                distance_data=self.distance_data,
                selected_genomes=self.selected_genomes,
                quality_scores=self.quality_scores,
                output_file=f.name,
                config=self.config
            )
            self.assertGreater(Path(f.name).stat().st_size, 1000)

    def test_create_visualization_size_deviation_column(self):
        """Test visualization with size_deviation_source column type."""
        config_with_size_dev = self.config.copy()
        config_with_size_dev['columns']['Assembly Size']['type'] = 'size_deviation_source'

        with tempfile.NamedTemporaryFile(suffix='.svg') as f:
            create_visualization(
                metadata=self.metadata,
                distance_data=self.distance_data,
                selected_genomes=self.selected_genomes,
                quality_scores=self.quality_scores,
                output_file=f.name,
                config=config_with_size_dev
            )
            self.assertGreater(Path(f.name).stat().st_size, 1000)

    def test_create_visualization_large_dataset(self):
        """Test visualization with larger dataset to check performance."""
        # Create larger dataset
        n_genomes = 20
        genome_ids = [f'genome{i:02d}' for i in range(n_genomes)]

        large_metadata = pd.DataFrame({
            'BUSCO [%S]': np.random.uniform(70, 98, n_genomes),
            'Assembly Size': np.random.uniform(3e6, 8e6, n_genomes),
            'Representative': np.random.choice([True, False], n_genomes),
            'Species': [f'Species_{i}' for i in range(n_genomes)]
        }, index=pd.Index(genome_ids, name='Identifier'))

        # Create distance matrix
        np.random.seed(42)
        distances = np.random.uniform(0, 0.3, (n_genomes, n_genomes))
        distances = (distances + distances.T) / 2  # Make symmetric
        np.fill_diagonal(distances, 0)

        large_distance = pd.DataFrame(distances, index=genome_ids, columns=genome_ids)

        large_selected = genome_ids[:5]  # Select first 5
        large_quality = {gid: np.random.uniform(0, 2) for gid in genome_ids}

        with tempfile.NamedTemporaryFile(suffix='.svg') as f:
            create_visualization(
                metadata=large_metadata,
                distance_data=large_distance,
                selected_genomes=large_selected,
                quality_scores=large_quality,
                output_file=f.name,
                config=self.config
            )
            self.assertGreater(Path(f.name).stat().st_size, 1000)


if __name__ == "__main__":
    unittest.main()
