"""
Tests for configuration validation and data loading edge cases.
"""


# Add src to path for imports
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd
from schema import SchemaError

sys.path.insert(0, str(Path(__file__).parent.parent))

from orthofinder_coreselector.data_loading import load_distance_matrix, load_metadata, validate_config


class TestConfigValidation(unittest.TestCase):
    """Test configuration validation functionality."""

    def test_validate_config_minimal(self):
        """Test validation of minimal valid config."""
        config = {
            'columns': {
                'test_col': {
                    'type': 'text'
                }
            }
        }

        validated = validate_config(config)

        self.assertIn('index', validated)
        self.assertEqual(validated['index'], 'Identifier')  # Default value
        self.assertIn('columns', validated)

    def test_validate_config_numeric_column(self):
        """Test validation of numeric column config."""
        config = {
            'columns': {
                'busco_score': {
                    'type': 'numeric',
                    'weight': 1.0,
                    'direction': 'higher_better',
                    'normalization': True,
                    'plot': True,
                    'display_unit': '%',
                    'scale_factor': 1.0,
                    'decimal_places': 2
                }
            }
        }

        validated = validate_config(config)
        col_config = validated['columns']['busco_score']

        self.assertEqual(col_config['type'], 'numeric')
        self.assertEqual(col_config['weight'], 1.0)
        self.assertEqual(col_config['direction'], 'higher_better')
        self.assertTrue(col_config['normalization'])
        self.assertTrue(col_config['plot'])

    def test_validate_config_boolean_column(self):
        """Test validation of boolean column config."""
        config = {
            'columns': {
                'representative': {
                    'type': 'boolean',
                    'true_value': True,
                    'weight': 0.5
                }
            }
        }

        validated = validate_config(config)
        col_config = validated['columns']['representative']

        self.assertEqual(col_config['type'], 'boolean')
        self.assertEqual(col_config['true_value'], True)
        self.assertEqual(col_config['weight'], 0.5)

    def test_validate_config_with_filter(self):
        """Test validation of config with filters."""
        config = {
            'columns': {
                'busco_score': {
                    'type': 'numeric',
                    'weight': 1.0,
                    'filter': {
                        'operator': '>=',
                        'value': 80.0
                    }
                }
            }
        }

        validated = validate_config(config)
        filter_config = validated['columns']['busco_score']['filter']

        self.assertEqual(filter_config['operator'], '>=')
        self.assertEqual(filter_config['value'], 80.0)

    def test_validate_config_missing_weight_numeric(self):
        """Test that numeric columns require weight."""
        config = {
            'columns': {
                'busco_score': {
                    'type': 'numeric'
                    # Missing weight
                }
            }
        }

        with self.assertRaises(SchemaError) as cm:
            validate_config(config)

        self.assertIn("must have 'weight'", str(cm.exception))

    def test_validate_config_missing_weight_boolean(self):
        """Test that boolean columns require weight."""
        config = {
            'columns': {
                'representative': {
                    'type': 'boolean',
                    'true_value': True
                    # Missing weight
                }
            }
        }

        with self.assertRaises(SchemaError) as cm:
            validate_config(config)

        self.assertIn("must have 'weight'", str(cm.exception))

    def test_validate_config_missing_true_value_boolean(self):
        """Test that boolean columns require true_value."""
        config = {
            'columns': {
                'representative': {
                    'type': 'boolean',
                    'weight': 1.0
                    # Missing true_value
                }
            }
        }

        with self.assertRaises(SchemaError) as cm:
            validate_config(config)

        self.assertIn("must specify 'true_value'", str(cm.exception))

    def test_validate_config_invalid_type(self):
        """Test validation with invalid column type."""
        config = {
            'columns': {
                'test_col': {
                    'type': 'invalid_type'
                }
            }
        }

        with self.assertRaises(SchemaError):
            validate_config(config)

    def test_validate_config_invalid_direction(self):
        """Test validation with invalid direction."""
        config = {
            'columns': {
                'busco_score': {
                    'type': 'numeric',
                    'weight': 1.0,
                    'direction': 'invalid_direction'
                }
            }
        }

        with self.assertRaises(SchemaError):
            validate_config(config)

    def test_validate_config_invalid_filter_operator(self):
        """Test validation with invalid filter operator."""
        config = {
            'columns': {
                'busco_score': {
                    'type': 'numeric',
                    'weight': 1.0,
                    'filter': {
                        'operator': 'invalid_op',
                        'value': 80.0
                    }
                }
            }
        }

        with self.assertRaises(SchemaError):
            validate_config(config)

    def test_validate_config_size_deviation_source(self):
        """Test validation of size_deviation_source column."""
        config = {
            'columns': {
                'cluster_id': {
                    'type': 'size_deviation_source',
                    'filter': {
                        'operator': '<=',
                        'value': 2.0
                    }
                }
            }
        }

        # Should not raise exception (size_deviation_source doesn't require weight)
        validated = validate_config(config)
        self.assertEqual(validated['columns']['cluster_id']['type'], 'size_deviation_source')


class TestDataLoadingEdgeCases(unittest.TestCase):
    """Test edge cases for data loading functions."""

    def test_load_metadata_empty_file(self):
        """Test loading empty metadata file."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write("Identifier\tSpecies\n")  # Header only
            temp_file = f.name

        try:
            with self.assertRaises(ValueError) as cm:
                load_metadata(temp_file, index_col="Identifier")

            self.assertIn("empty", str(cm.exception))
        finally:
            Path(temp_file).unlink()

    def test_load_metadata_missing_index_column(self):
        """Test loading metadata with missing index column."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write("Species\tBUSCO\n")
            f.write("E. coli\t95.0\n")
            temp_file = f.name

        try:
            with self.assertRaises(ValueError) as cm:
                load_metadata(temp_file, index_col="Identifier")

            self.assertIn("Error loading metadata file", str(cm.exception))
        finally:
            Path(temp_file).unlink()

    def test_load_metadata_malformed_file(self):
        """Test loading malformed metadata file."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write("This is not a valid TSV file\n")
            f.write("Random text without proper structure\n")
            temp_file = f.name

        try:
            with self.assertRaises(ValueError) as cm:
                load_metadata(temp_file, index_col="Identifier")

            self.assertIn("Error loading metadata", str(cm.exception))
        finally:
            Path(temp_file).unlink()

    def test_load_distance_matrix_csv_format(self):
        """Test loading distance matrix in CSV format."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write(",genome1,genome2,genome3\n")
            f.write("genome1,0.0,0.1,0.2\n")
            f.write("genome2,0.1,0.0,0.15\n")
            f.write("genome3,0.2,0.15,0.0\n")
            temp_file = f.name

        try:
            distance_data = load_distance_matrix(temp_file)

            self.assertIsInstance(distance_data, pd.DataFrame)
            self.assertEqual(len(distance_data), 3)
            self.assertEqual(len(distance_data.columns), 3)
            self.assertEqual(distance_data.loc['genome1', 'genome2'], 0.1)
        finally:
            Path(temp_file).unlink()

    def test_load_distance_matrix_tsv_format(self):
        """Test loading distance matrix in TSV format."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write("\tgenome1\tgenome2\tgenome3\n")
            f.write("genome1\t0.0\t0.1\t0.2\n")
            f.write("genome2\t0.1\t0.0\t0.15\n")
            f.write("genome3\t0.2\t0.15\t0.0\n")
            temp_file = f.name

        try:
            distance_data = load_distance_matrix(temp_file)

            self.assertIsInstance(distance_data, pd.DataFrame)
            self.assertEqual(len(distance_data), 3)
            self.assertEqual(len(distance_data.columns), 3)
            self.assertEqual(distance_data.loc['genome1', 'genome2'], 0.1)
        finally:
            Path(temp_file).unlink()

    def test_load_distance_matrix_malformed(self):
        """Test loading malformed distance matrix."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write("This is not a valid matrix file\n")
            temp_file = f.name

        try:
            # Malformed files may load as empty DataFrames rather than raising errors
            distance_data = load_distance_matrix(temp_file)
            # Should be empty or have issues
            self.assertTrue(len(distance_data) == 0 or len(distance_data.columns) == 0)
        finally:
            Path(temp_file).unlink()

    def test_load_metadata_successful(self):
        """Test successful metadata loading."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write("Identifier\tSpecies\tBUSCO [%S]\tRepresentative\n")
            f.write("genome1\tE. coli\t95.0\tTrue\n")
            f.write("genome2\tS. aureus\t90.0\tFalse\n")
            temp_file = f.name

        try:
            metadata = load_metadata(temp_file, index_col="Identifier")

            self.assertIsInstance(metadata, pd.DataFrame)
            self.assertEqual(len(metadata), 2)
            self.assertEqual(metadata.index.name, 'Identifier')
            self.assertIn('Species', metadata.columns)
            self.assertEqual(metadata.loc['genome1', 'Species'], 'E. coli')
        finally:
            Path(temp_file).unlink()


if __name__ == "__main__":
    unittest.main()
