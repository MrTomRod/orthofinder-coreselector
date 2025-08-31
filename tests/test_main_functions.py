"""
Tests for main.py functions.
"""

import logging

# Add src to path for imports
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))

from orthofinder_coreselector.data_loading import validate_config
from orthofinder_coreselector.main import (
    apply_filter,
    calculate_cluster_ids,
    calculate_final_scores,
    process_boolean,
    process_numeric,
    process_text,
    select_by_phylogenetic_clustering,
    select_core_set_api,
    validate_sample_consistency,
)


class TestMainFunctions(unittest.TestCase):
    """Test main.py functions."""

    def setUp(self):
        """Set up test data."""
        # Create test metadata
        self.metadata = pd.DataFrame({
            'BUSCO [%S]': [95.0, 90.0, 85.0, 80.0, 75.0],
            'Assembly Nr Scaffolds': [100, 200, 300, 400, 500],
            'Assembly Size': [5000000, 6000000, 7000000, 8000000, 9000000],
            'Representative': [True, False, True, False, False]
        }, index=pd.Index(['genome1', 'genome2', 'genome3', 'genome4', 'genome5'], name='Identifier'))

        # Create test distance matrix
        np.random.seed(42)  # For reproducible tests
        n = len(self.metadata)
        distances = np.random.rand(n, n) * 0.3
        distances = (distances + distances.T) / 2  # Make symmetric
        np.fill_diagonal(distances, 0)  # Zero diagonal

        self.distance_data = pd.DataFrame(
            distances,
            index=self.metadata.index,
            columns=self.metadata.index
        )

        # Test quality scores
        self.quality_scores = {
            'genome1': 1.5,
            'genome2': 0.8,
            'genome3': 1.2,
            'genome4': 0.5,
            'genome5': 0.3
        }

        # Test config
        self.config = {
            'index': 'Identifier',
            'columns': {
                'BUSCO [%S]': {
                    'type': 'numeric',
                    'weight': 1.0,
                    'direction': 'higher_better',
                    'plot': True
                },
                'Assembly Nr Scaffolds': {
                    'type': 'numeric',
                    'weight': 0.5,
                    'direction': 'lower_better',
                    'plot': True
                },
                'Representative': {
                    'type': 'boolean',
                    'true_value': True,
                    'weight': 1.0,
                    'plot': True
                }
            }
        }

        # Disable logging for tests
        logging.getLogger("orthofinder_coreselector").setLevel(logging.CRITICAL)

    def test_validate_sample_consistency_success(self):
        """Test successful sample consistency validation."""
        # Should not raise any exception
        validate_sample_consistency(self.metadata, self.distance_data)

    def test_validate_sample_consistency_mismatch(self):
        """Test sample consistency validation with mismatched samples."""
        # Create distance matrix with different samples
        bad_distance = self.distance_data.copy()
        bad_distance.index = ['genome1', 'genome2', 'genome3', 'genome4', 'genome_different']
        bad_distance.columns = bad_distance.index

        with self.assertRaises(ValueError) as cm:
            validate_sample_consistency(self.metadata, bad_distance)

        self.assertIn("Sample mismatch", str(cm.exception))

    def test_calculate_cluster_ids(self):
        """Test cluster ID calculation."""
        cluster_ids = calculate_cluster_ids(self.distance_data, target_core_size=3)

        self.assertIsInstance(cluster_ids, pd.Series)
        self.assertEqual(len(cluster_ids), len(self.distance_data))
        self.assertEqual(len(cluster_ids.unique()), 3)  # Should have 3 clusters
        self.assertTrue(all(0 <= cid < 3 for cid in cluster_ids))  # Cluster IDs 0-2

    def test_select_by_phylogenetic_clustering(self):
        """Test phylogenetic clustering selection."""
        selected, cluster_assignments = select_by_phylogenetic_clustering(
            self.quality_scores,
            self.distance_data,
            target_core_size=3,
            metadata=self.metadata
        )

        self.assertIsInstance(selected, list)
        self.assertIsInstance(cluster_assignments, dict)
        self.assertEqual(len(selected), 3)
        self.assertTrue(all(genome in self.metadata.index for genome in selected))

        # Check that cluster assignments cover all genomes
        self.assertEqual(len(cluster_assignments), len(self.metadata))

    def test_apply_filter_numeric(self):
        """Test applying numeric filters."""
        df = self.metadata.copy()
        df["__filter__"] = True

        # Test >= filter
        filter_config = {"operator": ">=", "value": 85.0}
        result = apply_filter(df, "BUSCO [%S]", filter_config)

        self.assertTrue((result[result["__filter__"]]["BUSCO [%S]"] >= 85.0).all())

        # Test < filter
        filter_config = {"operator": "<", "value": 85.0}
        result = apply_filter(df, "BUSCO [%S]", filter_config)

        filtered_genomes = result[result["__filter__"]]
        if len(filtered_genomes) > 0:
            self.assertTrue((filtered_genomes["BUSCO [%S]"] < 85.0).all())

    def test_apply_filter_boolean(self):
        """Test applying boolean filters."""
        df = self.metadata.copy()
        df["__filter__"] = True

        # Test == filter for boolean
        filter_config = {"operator": "==", "value": True}
        result = apply_filter(df, "Representative", filter_config)

        filtered_genomes = result[result["__filter__"]]
        if len(filtered_genomes) > 0:
            self.assertTrue(filtered_genomes["Representative"].all())

    def test_process_numeric(self):
        """Test numeric column processing."""
        df = self.metadata.copy()
        df["__filter__"] = True

        col_config = {
            'type': 'numeric',
            'weight': 1.0,
            'direction': 'higher_better',
            'normalization': True
        }

        result = process_numeric(df, "BUSCO [%S]", col_config)

        # Check that q-score column was created
        self.assertIn("BUSCO [%S]__q", result.columns)

        # Check that scores are z-scores (normalized with mean=0, std=1, then weighted)
        q_scores = result["BUSCO [%S]__q"]
        # Z-scores can be negative, check that they're reasonable values
        self.assertTrue(all(abs(q_scores) <= 10))  # Reasonable z-score range
        # Check that the scores have been weighted (weight=1.0 in config)
        self.assertAlmostEqual(q_scores.std(), 1.0, places=1)  # Should have std ~1 after z-score normalization

    def test_process_boolean(self):
        """Test boolean column processing."""
        df = self.metadata.copy()
        df["__filter__"] = True

        col_config = {
            'type': 'boolean',
            'true_value': True,
            'weight': 1.0
        }

        result = process_boolean(df, "Representative", col_config)

        # Check that q-score column was created
        self.assertIn("Representative__q", result.columns)

        # Check that True values get weight, False values get 0
        q_scores = result["Representative__q"]
        true_mask = result["Representative"] == True  # noqa: E712
        false_mask = result["Representative"] == False  # noqa: E712

        self.assertTrue((q_scores[true_mask] == 1.0).all())
        self.assertTrue((q_scores[false_mask] == 0.0).all())

    def test_process_text(self):
        """Test text column processing (should pass through unchanged)."""
        df = self.metadata.copy()
        df["text_col"] = ["A", "B", "C", "D", "E"]

        col_config = {'type': 'text'}
        result = process_text(df, "text_col", col_config)

        # Text processing should not add q-score columns
        self.assertNotIn("text_col__q", result.columns)
        # Original column should be unchanged
        pd.testing.assert_series_equal(result["text_col"], df["text_col"])

    def test_calculate_final_scores(self):
        """Test final score calculation."""
        df = pd.DataFrame({
            'col1__q': [0.8, 0.6, 0.4],
            'col2__q': [0.2, 0.3, 0.5],
            'other_col': [1, 2, 3]
        }, index=['genome1', 'genome2', 'genome3'])

        scores = calculate_final_scores(df)

        self.assertIsInstance(scores, dict)
        self.assertEqual(len(scores), 3)

        # Check that scores are sum of q-score columns
        expected_scores = {
            'genome1': 1.0,
            'genome2': 0.9,
            'genome3': 0.9
        }

        for genome, expected in expected_scores.items():
            self.assertAlmostEqual(scores[genome], expected, places=6)

    def test_calculate_final_scores_no_q_columns(self):
        """Test final score calculation with no q-score columns."""
        df = pd.DataFrame({
            'col1': [0.8, 0.6, 0.4],
            'col2': [0.2, 0.3, 0.5]
        }, index=['genome1', 'genome2', 'genome3'])

        with self.assertRaises(ValueError) as cm:
            calculate_final_scores(df)

        self.assertIn("No scoring columns found", str(cm.exception))

    def test_select_core_set_api(self):
        """Test the main API function."""
        config = validate_config(self.config)

        selected_genomes, quality_scores, cluster_assignments, filtered_genomes = select_core_set_api(
            metadata=self.metadata,
            distance_data=self.distance_data,
            target_core_size=3,
            config=config
        )

        self.assertIsInstance(selected_genomes, list)
        self.assertIsInstance(quality_scores, dict)
        self.assertIsInstance(cluster_assignments, dict)
        self.assertIsInstance(filtered_genomes, list)

        self.assertEqual(len(selected_genomes), 3)
        self.assertTrue(all(genome in self.metadata.index for genome in selected_genomes))
        self.assertEqual(len(quality_scores), len(filtered_genomes))
        self.assertEqual(len(cluster_assignments), len(self.metadata))

    def test_select_core_set_api_with_filters(self):
        """Test API function with filtering."""
        config_with_filter = self.config.copy()
        config_with_filter['columns']['BUSCO [%S]']['filter'] = {
            'operator': '>=',
            'value': 85.0
        }

        validated_config = validate_config(config_with_filter)

        selected_genomes, quality_scores, cluster_assignments, filtered_genomes = select_core_set_api(
            metadata=self.metadata,
            distance_data=self.distance_data,
            target_core_size=2,
            config=validated_config
        )

        # Should only have genomes with BUSCO >= 85%
        self.assertTrue(len(filtered_genomes) <= 3)  # genome1, genome2, genome3
        self.assertEqual(len(selected_genomes), 2)




if __name__ == "__main__":
    unittest.main()
