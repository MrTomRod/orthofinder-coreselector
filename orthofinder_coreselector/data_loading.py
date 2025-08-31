"""
Data loading utilities for OrthoFinder CoreSelector.
"""

from typing import Any, Dict

import pandas as pd
from schema import And, Or, Schema, SchemaError
from schema import Optional as SchemaOptional


def load_metadata(metadata_file: str, index_col: str) -> pd.DataFrame:
    """Load genome metadata from TSV file."""
    try:
        metadata = pd.read_csv(metadata_file, sep='\t')

        if len(metadata) == 0:
            raise ValueError("Metadata file is empty")

        assert index_col in metadata.columns, f"Index column '{index_col}' not found in metadata {metadata.columns=}"
        metadata = metadata.set_index(index_col)

        return metadata
    except Exception as e:
        raise ValueError(f"Error loading metadata file: {e}")


def load_distance_matrix(distance_matrix: str) -> pd.DataFrame:
    """Load ANI distance matrix."""
    try:
        # Try to detect file format and load accordingly
        if distance_matrix.endswith('.csv'):
            distance_data = pd.read_csv(distance_matrix, index_col=0)
        else:
            # Assume space/tab separated
            distance_data = pd.read_csv(distance_matrix, sep=None, engine='python', index_col=0)

        return distance_data
    except Exception as e:
        raise ValueError(f"Error loading distance matrix: {e}")


# Define the configuration schema
CONFIG_SCHEMA = Schema({
    SchemaOptional('index', default='Identifier'): str,
    'columns': {
        str: {  # Column name (any string)
            SchemaOptional('type', default='text'): And(str, lambda s: s in [
                'numeric', 'boolean', 'size_deviation_source', 'text'
            ]),
            SchemaOptional('weight'): And(Or(int, float), lambda n: isinstance(n, (int, float))),
            SchemaOptional('direction'): And(str, lambda s: s in ['higher_better', 'lower_better']),
            SchemaOptional('normalization', default=True): bool,
            SchemaOptional('true_value'): object,  # Any value for boolean columns
            SchemaOptional('plot', default=True): bool,
            SchemaOptional('display_unit', default=''): str,
            SchemaOptional('scale_factor', default=1.0): And(Or(int, float), lambda n: n > 0),
            SchemaOptional('decimal_places', default=2): And(int, lambda n: n >= 0),
            SchemaOptional('filter'): {
                'operator': And(str, lambda s: s in ['==', '!=', '>=', '<=', '>', '<', 'in', 'not_in']),
                'value': object
            }
        }
    }
})


def validate_config(config: Dict[str, Any]) -> Dict[str, Any]:
    """Validate and normalize JSON configuration against schema and metadata."""
    # Validate and normalize the config structure
    validated_config = CONFIG_SCHEMA.validate(config)

    # Validate that scoring columns have required fields
    for col_name, col_config in validated_config['columns'].items():
        col_type = col_config.get('type', 'text')

        if col_type in ['numeric', 'boolean']:
            if 'weight' not in col_config:
                raise SchemaError(f"Column '{col_name}' of type '{col_type}' must have 'weight'")

        # size_deviation_source can work without weight (filter-only mode)

        if col_type == 'boolean' and 'true_value' not in col_config:
            raise SchemaError(f"Boolean column '{col_name}' must specify 'true_value'")

    return validated_config
