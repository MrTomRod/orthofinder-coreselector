"""
Integration tests for OrthoFinder CoreSelector end-to-end workflows.
"""

import json
import logging

# Add src to path for imports
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))

from orthofinder_coreselector.data_loading import load_distance_matrix, load_metadata, validate_config
from orthofinder_coreselector.main import select_core_set_api
from orthofinder_coreselector.plots import create_visualization


class TestIntegration(unittest.TestCase):
    """Test complete workflows from start to finish."""

    def setUp(self):
        """Set up test data files."""
        # Disable logging for tests
        logging.getLogger("orthofinder_coreselector").setLevel(logging.CRITICAL)

        # Create realistic test data
        self.n_genomes = 10
        self.genome_ids = [f'GCA_{100000+i:06d}.1' for i in range(self.n_genomes)]

        # Create metadata
        np.random.seed(42)
        self.metadata_data = {
            'Identifier': self.genome_ids,
            'Species': [f'Species_{i//2}' for i in range(self.n_genomes)],  # Some duplicates
            'BUSCO [%S]': np.random.uniform(75, 98, self.n_genomes),
            'Assembly Size': np.random.uniform(3e6, 8e6, self.n_genomes),
            'Assembly Nr Scaffolds': np.random.randint(50, 500, self.n_genomes),
            'Representative': np.random.choice([True, False], self.n_genomes, p=[0.3, 0.7]),
            'N50': np.random.uniform(50000, 500000, self.n_genomes),
            'GC Content': np.random.uniform(35, 65, self.n_genomes)
        }

        # Create distance matrix
        distances = np.random.uniform(0, 0.4, (self.n_genomes, self.n_genomes))
        distances = (distances + distances.T) / 2  # Make symmetric
        np.fill_diagonal(distances, 0)  # Zero diagonal
        self.distance_matrix = distances

        # Create comprehensive config
        self.config = {
            'index': 'Identifier',
            'columns': {
                'BUSCO [%S]': {
                    'type': 'numeric',
                    'weight': 2.0,
                    'direction': 'higher_better',
                    'plot': True,
                    'display_unit': '%',
                    'decimal_places': 1,
                    'filter': {
                        'operator': '>=',
                        'value': 80.0
                    }
                },
                'Assembly Size': {
                    'type': 'numeric',
                    'weight': 1.0,
                    'direction': 'higher_better',
                    'plot': True,
                    'display_unit': 'Mbp',
                    'scale_factor': 1000000,
                    'decimal_places': 1
                },
                'Assembly Nr Scaffolds': {
                    'type': 'numeric',
                    'weight': 0.5,
                    'direction': 'lower_better',
                    'plot': True,
                    'display_unit': 'scaffolds',
                    'decimal_places': 0
                },
                'Representative': {
                    'type': 'boolean',
                    'true_value': True,
                    'weight': 1.5,
                    'plot': True
                },
                'N50': {
                    'type': 'numeric',
                    'weight': 1.0,
                    'direction': 'higher_better',
                    'plot': True,
                    'display_unit': 'bp',
                    'decimal_places': 0
                },
                'Species': {
                    'type': 'text',
                    'plot': True
                },
                'GC Content': {
                    'type': 'numeric',
                    'weight': 0.2,
                    'direction': 'higher_better',
                    'plot': False  # Don't plot this one
                }
            }
        }

    def create_temp_files(self):
        """Create temporary files for testing."""
        # Create metadata file
        metadata_df = pd.DataFrame(self.metadata_data)
        metadata_file = tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False)
        metadata_df.to_csv(metadata_file.name, sep='\t', index=False)
        metadata_file.close()

        # Create distance matrix file
        distance_df = pd.DataFrame(
            self.distance_matrix,
            index=self.genome_ids,
            columns=self.genome_ids
        )
        distance_file = tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False)
        distance_df.to_csv(distance_file.name)
        distance_file.close()

        # Create config file
        config_file = tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False)
        json.dump(self.config, config_file, indent=2)
        config_file.close()

        return metadata_file.name, distance_file.name, config_file.name

    def cleanup_temp_files(self, *file_paths):
        """Clean up temporary files."""
        for file_path in file_paths:
            if file_path and Path(file_path).exists():
                Path(file_path).unlink()

    def test_complete_workflow_small_target(self):
        """Test complete workflow with small target core size."""
        metadata_file, distance_file, config_file = self.create_temp_files()

        try:
            # Load data
            metadata = load_metadata(metadata_file, index_col='Identifier')
            distance_data = load_distance_matrix(distance_file)

            with open(config_file) as f:
                config = json.load(f)
            validated_config = validate_config(config)

            # Run core selection
            selected_genomes, quality_scores, cluster_assignments, filtered_genomes = select_core_set_api(
                metadata=metadata,
                distance_data=distance_data,
                target_core_size=3,
                config=validated_config
            )

            # Validate results
            self.assertIsInstance(selected_genomes, list)
            self.assertEqual(len(selected_genomes), 3)
            self.assertTrue(all(genome in metadata.index for genome in selected_genomes))

            # Check that filtering worked (should have some genomes with BUSCO >= 80%)
            self.assertGreater(len(filtered_genomes), 0)
            for genome in filtered_genomes:
                self.assertGreaterEqual(metadata.loc[genome, 'BUSCO [%S]'], 80.0)

            # Test visualization
            with tempfile.NamedTemporaryFile(suffix='.svg') as viz_file:
                create_visualization(
                    metadata=metadata,
                    distance_data=distance_data,
                    selected_genomes=selected_genomes,
                    quality_scores=quality_scores,
                    output_file=viz_file.name,
                    config=validated_config,
                    cluster_assignments=cluster_assignments,
                    filtered_genomes=filtered_genomes
                )

                # Check that visualization file has content (not empty)
                self.assertGreater(Path(viz_file.name).stat().st_size, 1000)

        finally:
            self.cleanup_temp_files(metadata_file, distance_file, config_file)

    def test_complete_workflow_large_target(self):
        """Test complete workflow with target larger than available genomes."""
        metadata_file, distance_file, config_file = self.create_temp_files()

        try:
            # Load data
            metadata = load_metadata(metadata_file, index_col='Identifier')
            distance_data = load_distance_matrix(distance_file)

            with open(config_file) as f:
                config = json.load(f)
            validated_config = validate_config(config)

            # Request more genomes than available after filtering - should raise ValueError
            with self.assertRaises(ValueError) as cm:
                select_core_set_api(
                    metadata=metadata,
                    distance_data=distance_data,
                    target_core_size=20,  # More than total genomes
                    config=validated_config
                )

            # Check error message is clear
            self.assertIn("Only", str(cm.exception))
            self.assertIn("available after filtering", str(cm.exception))

        finally:
            self.cleanup_temp_files(metadata_file, distance_file, config_file)



    def test_workflow_boolean_only_scoring(self):
        """Test workflow using only boolean columns for scoring."""
        boolean_config = {
            'index': 'Identifier',
            'columns': {
                'Representative': {
                    'type': 'boolean',
                    'true_value': True,
                    'weight': 1.0,
                    'plot': True
                },
                'Species': {
                    'type': 'text',
                    'plot': True
                },
                'BUSCO [%S]': {
                    'type': 'text',  # Use text type to avoid weight requirement
                    'plot': True  # Plot but don't score
                }
            }
        }

        metadata_file, distance_file, _ = self.create_temp_files()
        config_file = tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False)
        json.dump(boolean_config, config_file, indent=2)
        config_file.close()

        try:
            # Load data
            metadata = load_metadata(metadata_file, index_col='Identifier')
            distance_data = load_distance_matrix(distance_file)
            validated_config = validate_config(boolean_config)

            # Run core selection
            selected_genomes, quality_scores, cluster_assignments, filtered_genomes = select_core_set_api(
                metadata=metadata,
                distance_data=distance_data,
                target_core_size=3,
                config=validated_config
            )

            # Should work with only boolean scoring
            self.assertEqual(len(selected_genomes), 3)

            # Quality scores should only reflect Representative column
            for genome, score in quality_scores.items():
                expected_score = 1.0 if metadata.loc[genome, 'Representative'] else 0.0
                self.assertEqual(score, expected_score)

        finally:
            self.cleanup_temp_files(metadata_file, distance_file, config_file.name)

    def test_workflow_with_missing_data(self):
        """Test workflow robustness with missing data."""
        # Create metadata with some missing values
        metadata_with_missing = self.metadata_data.copy()
        # Introduce some NaN values
        for i in [1, 3, 5]:
            metadata_with_missing['BUSCO [%S]'][i] = None
            metadata_with_missing['N50'][i] = None

        metadata_file = tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False)
        pd.DataFrame(metadata_with_missing).to_csv(metadata_file.name, sep='\t', index=False)
        metadata_file.close()

        _, distance_file, config_file = self.create_temp_files()

        try:
            # Load data
            metadata = load_metadata(metadata_file.name, index_col='Identifier')
            distance_data = load_distance_matrix(distance_file)

            with open(config_file) as f:
                config = json.load(f)
            validated_config = validate_config(config)

            # Run core selection - should handle missing data gracefully
            selected_genomes, quality_scores, cluster_assignments, filtered_genomes = select_core_set_api(
                metadata=metadata,
                distance_data=distance_data,
                target_core_size=3,
                config=validated_config
            )

            # Should still work
            self.assertIsInstance(selected_genomes, list)
            self.assertGreater(len(selected_genomes), 0)

        finally:
            self.cleanup_temp_files(metadata_file.name, distance_file, config_file)

    def test_reproducibility(self):
        """Test that results are reproducible with same input."""
        metadata_file, distance_file, config_file = self.create_temp_files()

        try:
            # Load data
            metadata = load_metadata(metadata_file, index_col='Identifier')
            distance_data = load_distance_matrix(distance_file)

            with open(config_file) as f:
                config = json.load(f)
            validated_config = validate_config(config)

            # Run twice with same parameters
            result1 = select_core_set_api(
                metadata=metadata,
                distance_data=distance_data,
                target_core_size=4,
                config=validated_config
            )

            result2 = select_core_set_api(
                metadata=metadata,
                distance_data=distance_data,
                target_core_size=4,
                config=validated_config
            )

            # Results should be identical
            self.assertEqual(result1[0], result2[0])  # selected_genomes
            self.assertEqual(result1[1], result2[1])  # quality_scores
            self.assertEqual(result1[2], result2[2])  # cluster_assignments
            self.assertEqual(result1[3], result2[3])  # filtered_genomes

        finally:
            self.cleanup_temp_files(metadata_file, distance_file, config_file)


if __name__ == "__main__":
    unittest.main()
