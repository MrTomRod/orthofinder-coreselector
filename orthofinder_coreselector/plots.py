"""
Visualization utilities for OrthoFinder CoreSelector.
"""

import logging
from typing import Dict, List

import matplotlib
import pandas as pd

matplotlib.use('SVG')  # Use SVG backend
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform

# Configure matplotlib for SVG text rendering
plt.rcParams['svg.fonttype'] = 'none'  # Ensure text is saved as text, not paths

# Set up module-level logger
logger = logging.getLogger(__name__)


def create_visualization(
    metadata: pd.DataFrame,
    distance_data: pd.DataFrame,
    selected_genomes: List[str],
    quality_scores: Dict[str, float],
    output_file: str,
    cluster_assignments: Dict[str, int] = None,
    config: Dict = None,
    filtered_genomes: List[str] = None
):
    """Create SVG visualization with tree and metadata plots to the right."""
    # Get common genomes between metadata and distance matrix
    common_genomes = list(set(metadata.index.tolist()) &
                        set(distance_data.index.tolist()))

    if len(common_genomes) == 0:
        raise ValueError("No common genomes found between metadata and distance matrix")

    # Ensure all selected genomes are included in visualization
    selected_in_common = [g for g in selected_genomes if g in common_genomes]
    logger.info(f"   Selected genomes in common set: {len(selected_in_common)}/{len(selected_genomes)}")

    # Use all common genomes for visualization
    logger.info(f"   Creating visualization with {len(common_genomes)} genomes")

    # Use distance matrix directly for clustering
    subset_distance = distance_data.loc[common_genomes, common_genomes]

    # Perform hierarchical clustering
    condensed_matrix = squareform(subset_distance.values)
    linkage_matrix = hierarchy.linkage(condensed_matrix, method='average')

    genome_names = common_genomes

    # Get the order of genomes from clustering
    # We'll get the actual leaf positions from the real dendrogram
    ordered_genomes = genome_names

    # Calculate dynamic height based on number of leaves
    n_leaves = len(ordered_genomes)
    tree_height = max(8, n_leaves * 0.3)  # Minimum 8 inches, 0.3 inches per leaf

    # Set up the figure with dynamic height
    fig_width = 20
    fig_height = tree_height
    fig = plt.figure(figsize=(fig_width, fig_height))

    # Create grid for tree and metadata plots
    # Tree takes 40% of width, metadata plots take 60%
    # Add more space between subplots for leaf labels
    gs = fig.add_gridspec(1, 2, width_ratios=[0.4, 0.6], wspace=0.3)

    # Tree subplot
    ax_tree = fig.add_subplot(gs[0])

    # Create dendrogram with black branches
    dendro_data = hierarchy.dendrogram(
        linkage_matrix,
        labels=ordered_genomes,
        ax=ax_tree,
        orientation='left',
        leaf_font_size=8,
        distance_sort='descending',
        link_color_func=lambda x: 'black'
    )

    # Color leaf labels based on cluster assignments if provided
    if cluster_assignments:
        import matplotlib.colors as mcolors
        colors = list(mcolors.TABLEAU_COLORS.values())  # Use Tableau colors for good contrast

        # Get the leaf labels (text objects)
        leaf_labels = ax_tree.get_yticklabels()

        # Color each leaf label based on its cluster assignment
        for i, label in enumerate(leaf_labels):
            genome = label.get_text()
            if genome in cluster_assignments:
                cluster_id = cluster_assignments[genome]
                color = colors[cluster_id % len(colors)]
                label.set_color(color)

                # Font weight based on filter status
                if filtered_genomes and genome in filtered_genomes:
                    label.set_fontweight('bold')  # Non-filtered genomes are bold
                else:
                    label.set_fontweight('normal')  # Filtered out genomes are regular
    ax_tree.set_title('Hierarchical Clustering Tree (ANI-based)', fontsize=14, fontweight='bold')
    ax_tree.set_xlabel('Distance', fontsize=12)

    # Get the actual leaf order from the dendrogram
    leaf_order = dendro_data['leaves']

    # The dendrogram reorders the genomes - get the new order
    # Ensure we don't go out of bounds
    ordered_genomes_dendro = []
    for i in leaf_order:
        if i < len(ordered_genomes):
            ordered_genomes_dendro.append(ordered_genomes[i])
        else:
            raise ValueError(f"Leaf index {i} out of bounds for {len(ordered_genomes)} genomes. "
            "This indicates a bug in the clustering logic.")

    # In the dendrogram plot, leaves appear at y = 0, 1, 2, ..., n-1 in the reordered sequence
    # So we need to find where each original genome appears in this new order
    actual_y_positions = []
    for genome in ordered_genomes:
        # Find where this genome appears in the dendrogram order
        try:
            dendro_position = ordered_genomes_dendro.index(genome)
            actual_y_positions.append(dendro_position)
        except ValueError:
            raise ValueError(f"Genome {genome} not found in dendrogram order. This indicates a data consistency bug.")

    # Highlight selected genomes in the tree using actual positions
    for i, genome in enumerate(ordered_genomes):
        if genome in selected_genomes:
            # Use the actual y-position where this genome appears in the dendrogram
            y_pos = actual_y_positions[i]
            ax_tree.axhline(y=5+y_pos*10, color='red', alpha=0.3, linewidth=2)

    # Create metadata plots to the right using the same positions (only if we have columns to plot)
    if 'columns' in config and config['columns']:
        create_metadata_plots_right(
            fig=fig,
            gs_metadata=gs[1],
            metadata=metadata,
            ordered_genomes=ordered_genomes,
            selected_genomes=selected_genomes,
            quality_scores=quality_scores,
            actual_y_positions=actual_y_positions,
            config=config,
            cluster_assignments=cluster_assignments,
            filtered_genomes=filtered_genomes
        )

    plt.savefig(output_file, format='svg', bbox_inches='tight', metadata={'Creator': 'OrthoFinder CoreSelector'})
    plt.close()

    logger.info(f"   Visualization saved to: {output_file}")




def create_metadata_plots_right(
    fig,
    gs_metadata,
    metadata: pd.DataFrame,
    ordered_genomes: List[str],
    selected_genomes: List[str],
    quality_scores: Dict[str, float],
    actual_y_positions: List[int],
    config: Dict = None,
    cluster_assignments: Dict[str, int] = None,
    filtered_genomes: List[str] = None

):
    """Create metadata plots to the right of the tree, aligned with leaf positions."""
    # Define colors for cluster assignments
    import matplotlib.colors as mcolors
    colors = list(mcolors.TABLEAU_COLORS.values())  # Use Tableau colors for good contrast

    # Get columns to plot from config
    metadata_columns = []
    for col_name, col_config in config['columns'].items():
        # Only plot columns that have plot=true and exist in metadata
        if col_config.get('plot', False) and col_name in metadata.columns:
            col_type = col_config.get('type', 'text')

            # Create appropriate transform function and title based on column type
            if col_type == 'numeric':
                # Use config metadata for units/formatting if available
                display_unit = col_config.get('display_unit', '')
                if display_unit:
                    title = f"{col_name} ({display_unit})"
                else:
                    title = col_name

                # Use scale factor from config if available
                scale_factor = col_config.get('scale_factor', 1.0)
                # Create lambda with proper closure
                transform_func = (lambda sf: lambda x: float(x) / sf if pd.notna(x) and x != 'None' else 0)(scale_factor)

            elif col_type == 'size_deviation_source':
                # For size deviation, assume it's assembly size and show in MB
                display_unit = col_config.get('display_unit', 'MB')
                title = f"{col_name} ({display_unit})"
                scale_factor = col_config.get('scale_factor', 1e6)  # Default to MB scale
                # Create lambda with proper closure
                transform_func = (lambda sf: lambda x: float(x) / sf if pd.notna(x) and x != 'None' else 0)(scale_factor)
            elif col_type == 'boolean':
                title = col_name
                true_value = col_config.get('true_value', True)
                # Create lambda with proper closure
                transform_func = (lambda tv: lambda x: 1 if str(x) == str(tv) else 0)(true_value)
            elif col_type == 'text':
                title = col_name
                transform_func = lambda x: str(x) if pd.notna(x) else ""
            else:
                title = col_name
                transform_func = lambda x: x

            metadata_columns.append((col_name, title, transform_func, col_type))

    n_plots = len(metadata_columns)
    n_leaves = len(ordered_genomes)

    if n_plots == 0:
        raise ValueError("No metadata columns marked for plotting")
    if n_leaves < 3:
        raise ValueError("Less than 3 genomes in the dataset")

    # Create subplots for each metadata column
    metadata_axes = gs_metadata.subgridspec(1, n_plots, wspace=0.05)

    for plot_idx, (col_name, plot_title, transform_func, col_type) in enumerate(metadata_columns):
        ax = fig.add_subplot(metadata_axes[plot_idx])

        # Handle text columns differently from numeric/boolean columns
        if col_type == 'text':
            # For text columns, just show the text values as labels
            ax.set_xlim(0, 1)
            for i, genome in enumerate(ordered_genomes):
                if genome in metadata.index:
                    row = metadata.loc[genome]
                    text_value = transform_func(row[col_name])

                    # Color based on cluster assignment
                    if cluster_assignments and genome in cluster_assignments:
                        cluster_id = cluster_assignments[genome]
                        color = colors[cluster_id % len(colors)]
                    else:
                        color = 'black'

                    # Text columns use normal font weight (no bold/regular variation)
                    ax.text(0.5, actual_y_positions[i], text_value,
                           ha='center', va='center', fontsize=8,
                           color=color, fontweight='normal')
                else:
                    ax.text(0.5, actual_y_positions[i], "",
                           ha='center', va='center', fontsize=8, color='gray')

            # Remove x-axis for text columns since there are no bars
            ax.set_xticks([])
            ax.spines['bottom'].set_visible(False)

        else:
            # Prepare data for numeric/boolean columns
            values = []
            colors = []

            for genome in ordered_genomes:
                if genome in metadata.index:
                    row = metadata.loc[genome]
                    value = transform_func(row[col_name])
                    values.append(value)

                    # Color based on selection status and column type
                    if col_type == 'boolean':
                        # For boolean columns, use green for True, red for False
                        colors.append('green' if value == 1 else 'red')
                        values[-1] = 1  # Make all bars the same height, we only care about color
                    else:
                        # For other columns, use red for selected, lightblue for not selected
                        if genome in selected_genomes:
                            colors.append('red')
                        else:
                            colors.append('lightblue')
                else:
                    values.append(0)
                    colors.append('lightgray')

            # Create horizontal bar plot (aligned with tree leaves)
            # Use the actual y-positions from the dendrogram
            ax.barh(actual_y_positions, values, color=colors, alpha=0.7, height=0.8)

        # Customize the plot
        ax.set_title(plot_title, fontsize=12, fontweight='bold')
        ax.set_xlabel('Value', fontsize=10)

        # Set y-axis to match tree leaf positions
        ax.set_ylim(min(actual_y_positions) - 0.5, max(actual_y_positions) + 0.5)
        ax.set_yticks(actual_y_positions)
        ax.set_yticklabels([])  # Hide y-axis labels to avoid clutter

        # Add grid
        ax.grid(True, alpha=0.3, axis='x')

        # Add value annotations for non-text columns (text columns already have their labels)
        if col_type != 'text':
            for i, genome in enumerate(ordered_genomes):
                if genome in metadata.index:
                    # For boolean plots, show quality score
                    if col_type == 'boolean':
                        if genome in quality_scores:
                            score = quality_scores[genome]
                            ax.text(0.5, actual_y_positions[i], f'Q:{score:.2f}',
                                   ha='center', va='center', fontsize=8)
                    else:
                        # For other plots, show the actual value
                        value = values[i]
                        # Format the value based on config settings
                        # Get current column config from the metadata_columns tuple
                        current_col_config = config['columns'].get(col_name, {})
                        display_unit = current_col_config.get('display_unit', '')
                        decimal_places = current_col_config.get('decimal_places', 2)

                        if col_type == 'numeric' or col_type == 'size_deviation_source':
                            if isinstance(value, float):
                                if value == int(value):  # Integer values
                                    text_value = f'{int(value)}'
                                else:  # Float values
                                    text_value = f'{value:.{decimal_places}f}'
                            else:
                                text_value = str(value)

                            # Add unit if specified
                            if display_unit:
                                text_value += display_unit
                        else:
                            text_value = str(value)

                        # Position text at the end of the bar
                        text_x = max(values) * 0.5 if max(values) > 0 else 0.5
                        ax.text(text_x, actual_y_positions[i], text_value,
                               ha='center', va='center',
                               fontsize=7, color='black')

        # Remove top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Only show x-axis labels for the leftmost plot
        if plot_idx == 0:
            ax.set_ylabel('Genomes', fontsize=10)
        else:
            ax.set_ylabel('')
