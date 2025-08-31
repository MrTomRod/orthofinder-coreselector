"""
Tests for visualization functionality.
"""

# Add src to path for imports
import sys
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))

from orthofinder_coreselector.plots import create_visualization


class TestVisualization(unittest.TestCase):
    """Test visualization functionality."""

    def setUp(self):
        """Set up test data."""
        # Create simple test metadata with matching genome IDs
        self.metadata = pd.DataFrame({
            'BUSCO [%S]': [95.0, 90.0, 85.0],
            'Assembly Nr Scaffolds': [100, 200, 300],
            'Assembly Size': [5000000, 6000000, 7000000],
            'Representative': [True, False, False]
        }, index=pd.Index(['genome1', 'genome2', 'genome3'], name='Identifier'))

        # Create simple distance matrix with matching genome IDs
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

        # Minimal config for testing
        self.config = {
            'columns': {
                'BUSCO [%S]': {
                    'type': 'numeric',
                    'weight': 1.0,
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

    def test_visualization_creation(self):
        """Test visualization creation."""
        output_file = Path(__file__).parent / "test_output.svg"

        try:
            create_visualization(
                metadata=self.metadata,
                distance_data=self.distance_data,
                selected_genomes=self.selected_genomes,
                quality_scores=self.quality_scores,
                output_file=str(output_file),
                config=self.config
            )

            # Check that output file was created
            self.assertTrue(output_file.exists())

        finally:
            # Clean up
            if output_file.exists():
                output_file.unlink()


if __name__ == "__main__":
    unittest.main()
