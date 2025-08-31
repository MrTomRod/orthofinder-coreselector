"""
OrthoFinder CoreSelector: Automated Core Set Selector for OrthoFinder v3

This tool selects an optimal core set of genomes for OrthoFinder analysis using
a flexible, JSON-configured workflow:
1. Loading and validating configuration, metadata, and distance matrices
2. Processing columns according to configuration (scoring, filtering, plotting)
3. Selecting genomes using phylogenetic clustering for maximum diversity
4. Generating comprehensive visualizations

The system supports arbitrary column types (numeric, boolean, size_deviation_source, text)
with configurable weights, filters, and visualization options.
"""

__version__ = "0.1.0"

from .data_loading import load_distance_matrix, load_metadata, validate_config
from .init_orthofinder_arx import init_orthofinder_arx
from .main import select_core_set, select_core_set_api
from .plots import create_visualization

__all__ = [
    "load_metadata",
    "load_distance_matrix",
    "validate_config",
    "create_visualization",
    "select_core_set",
    "select_core_set_api",
    "init_orthofinder_arx"
]
