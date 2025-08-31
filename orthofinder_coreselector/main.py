"""
OrthoFinder CoreSelector: Automated Core Set Selector for OrthoFinder v3

This tool selects an optimal core set of genomes for OrthoFinder analysis using
a flexible, JSON-configured workflow:
1. Loading and validating configuration, metadata, and distance matrices
2. Processing columns according to configuration (scoring, filtering, plotting)
3. Selecting genomes using phylogenetic clustering for maximum diversity
4. Generating comprehensive visualizations

Usage:
    # Minimal mode (random selection within phylogenetic clusters)
    orthofinder-coreselector <distance_matrix> <target_core_size> <output_file> <visualization_file>
                            [--log-level INFO] [--max-genomes N]

    # Quality mode (metadata-based selection within clusters)
    orthofinder-coreselector <distance_matrix> <target_core_size> <output_file> <visualization_file>
                            --metadata-file <metadata_file> --config-file <config_file>
                            [--log-level INFO] [--max-genomes N]
"""

import json
import logging
import sys
from typing import Any, Dict, List, Optional, Tuple

import fire
import numpy as np
import pandas as pd
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform

from .data_loading import load_distance_matrix, load_metadata, validate_config
from .plots import create_visualization

# Set up module-level logger
logger = logging.getLogger(__name__)


def validate_sample_consistency(
    metadata: pd.DataFrame,
    distance_data: pd.DataFrame,
) -> None:
    """
    Validate that metadata and distance matrix contain the same samples.

    Args:
        metadata: Metadata DataFrame
        distance_data: Distance matrix DataFrame

    Raises:
        ValueError: If samples don't match between files
    """
    metadata_genomes = set(metadata.index.tolist())
    distance_genomes = set(distance_data.index.tolist())

    # Check if all samples are present in both files
    if metadata_genomes != distance_genomes:
        missing_in_distance = metadata_genomes - distance_genomes
        missing_in_metadata = distance_genomes - metadata_genomes

        error_msg = "❌ Sample mismatch between metadata and distance matrix files:\n"

        if missing_in_distance:
            error_msg += f"   {len(missing_in_distance)} samples in metadata but missing from distance matrix:\n"
            for sample in sorted(list(missing_in_distance))[:10]:  # Show first 10
                error_msg += f"     - {sample}\n"
            if len(missing_in_distance) > 10:
                error_msg += f"     ... and {len(missing_in_distance) - 10} more\n"

        if missing_in_metadata:
            error_msg += f"   {len(missing_in_metadata)} samples in distance matrix but missing from metadata:\n"
            for sample in sorted(list(missing_in_metadata))[:10]:  # Show first 10
                error_msg += f"     - {sample}\n"
            if len(missing_in_metadata) > 10:
                error_msg += f"     ... and {len(missing_in_metadata) - 10} more\n"

        error_msg += f"\n   Total samples in metadata: {len(metadata_genomes)}"
        error_msg += f"\n   Total samples in distance matrix: {len(distance_genomes)}"
        error_msg += f"\n   Common samples: {len(metadata_genomes & distance_genomes)}"
        error_msg += "\n\nPlease ensure both files contain the same samples."

        raise ValueError(error_msg)

    logger.info(f"✅ All {len(metadata_genomes)} samples are present in both metadata and distance matrix files")


def select_by_phylogenetic_clustering(
    quality_scores: Dict[str, float],
    distance_data: pd.DataFrame,
    target_core_size: int,
    metadata: pd.DataFrame,
) -> tuple[List[str], Dict[str, int]]:
    """
    Select genomes by creating phylogenetic clusters and picking best from each.

    Args:
        quality_scores: Dictionary mapping genome IDs to quality scores
        distance_data: Distance matrix DataFrame
        target_core_size: Number of genomes to select
        metadata: Metadata DataFrame for additional genome information

    Returns:
        Tuple of (selected genome IDs, cluster assignments dictionary)
    """
    logger.info(f"   Creating {target_core_size} phylogenetic clusters...")

    # Handle edge cases: raise error if there are no genomes or only one genome
    if len(distance_data) < 3:
        raise ValueError("Less than 3 genomes in the dataset")

    # Create hierarchical clustering
    condensed_matrix = squareform(distance_data.values)
    linkage_matrix = hierarchy.linkage(condensed_matrix, method='average')

    # Cut tree into target_core_size clusters
    cluster_labels = hierarchy.cut_tree(linkage_matrix, n_clusters=target_core_size)

    # Group genomes by cluster
    clusters = {}
    for i, genome in enumerate(distance_data.index):
        cluster_id = int(cluster_labels.flat[i])  # Use proper numpy indexing
        if cluster_id not in clusters:
            clusters[cluster_id] = []
        clusters[cluster_id].append(genome)

    logger.info(f"   Created {len(clusters)} clusters")
    cluster_sizes = [len(genomes) for genomes in clusters.values()]
    logger.info(f"   Cluster sizes: min={min(cluster_sizes)}, max={max(cluster_sizes)}, avg={np.mean(cluster_sizes):.1f}")

    # Select best genome from each cluster
    selected = []
    for cluster_id, genomes in clusters.items():
        # Find the genome with highest quality score in this cluster
        best_genome = max(genomes, key=lambda g: quality_scores.get(g, 0))
        selected.append(best_genome)

        best_score = quality_scores.get(best_genome, 0)
        logger.debug(f"   Cluster {cluster_id + 1}/{target_core_size}: selected {best_genome} (Q={best_score:.3f}) from {len(genomes)} genomes")
    logger.debug("")  # Spacing for readability

    # Create cluster assignments dictionary
    cluster_assignments = {}
    for i, genome in enumerate(distance_data.index):
        # Use proper numpy indexing to get scalar value
        cluster_id = int(cluster_labels.flat[i])
        cluster_assignments[genome] = cluster_id

    return selected, cluster_assignments



def calculate_cluster_ids(distance_data: pd.DataFrame, target_core_size: int) -> pd.Series:
    """Calculate cluster IDs for all genomes."""
    condensed_matrix = squareform(distance_data.values)
    linkage_matrix = hierarchy.linkage(condensed_matrix, method='average')
    cluster_labels = hierarchy.cut_tree(linkage_matrix, n_clusters=target_core_size)

    return pd.Series(
        cluster_labels.flat,
        index=distance_data.index,
        name="cluster_id"
    )


def apply_filter(df: pd.DataFrame, col_name: str, filter_config: Dict[str, Any]) -> pd.DataFrame:
    """Apply filter to column and update __filter__ column."""
    operator = filter_config["operator"]
    value = filter_config["value"]

    # Count genomes before filtering
    before_count = df["__filter__"].sum()

    if operator == "==":
        mask = df[col_name] == value
    elif operator == "!=":
        mask = df[col_name] != value
    elif operator == ">=":
        mask = df[col_name].astype(float) >= value
    elif operator == "<=":
        mask = df[col_name].astype(float) <= value
    elif operator == ">":
        mask = df[col_name].astype(float) > value
    elif operator == "<":
        mask = df[col_name].astype(float) < value
    elif operator == "in":
        mask = df[col_name].isin(value)
    elif operator == "not_in":
        mask = ~df[col_name].isin(value)
    else:
        raise ValueError(f"Unknown operator: {operator}")

    df["__filter__"] = df["__filter__"] & mask

    # Count genomes after filtering and report
    after_count = df["__filter__"].sum()
    removed_count = before_count - after_count

    logger.info(f"     Filter {col_name} {operator} {value}: removed {removed_count} genomes ({before_count} → {after_count})")

    return df


def process_numeric(df: pd.DataFrame, col_name: str, col_config: Dict[str, Any]) -> pd.DataFrame:
    """Process numeric column: normalize, score, filter."""
    # Handle missing values and convert to float
    values = df[col_name].replace('None', np.nan).astype(float)
    valid_mask = ~values.isna()

    if valid_mask.sum() == 0:
        raise ValueError(f"No valid numeric values in column '{col_name}'")

    # Check if normalization should be applied
    use_normalization = col_config.get("normalization", True)

    if use_normalization:
        # Calculate z-scores for valid values (existing behavior)
        valid_values = values[valid_mask]
        mean_val = valid_values.mean()
        std_val = valid_values.std()

        if std_val == 0:
            df[f"{col_name}__q"] = 0.0
        else:
            z_scores = (values - mean_val) / std_val

            # Apply direction
            direction = col_config.get("direction", "higher_better")
            if direction == "lower_better":
                z_scores = -z_scores

            # Apply weight
            weight = col_config.get("weight", 0.0)
            df[f"{col_name}__q"] = z_scores * weight

            # Set NaN scores to 0
            df[f"{col_name}__q"] = df[f"{col_name}__q"].fillna(0.0)
    else:
        # Direct scoring without normalization (new behavior)
        scores = values.copy()

        # Apply direction
        direction = col_config.get("direction", "higher_better")
        if direction == "lower_better":
            # For lower_better, invert the scores (1 - score for 0-1 range)
            max_val = scores.max()
            min_val = scores.min()
            if max_val <= 1.0 and min_val >= 0.0:
                # Assume 0-1 range, invert it
                scores = 1.0 - scores
            else:
                # For other ranges, use negative values
                scores = -scores

        # Apply weight
        weight = col_config.get("weight", 0.0)
        df[f"{col_name}__q"] = scores * weight

        # Set NaN scores to 0
        df[f"{col_name}__q"] = df[f"{col_name}__q"].fillna(0.0)

    # Apply filter if specified
    if "filter" in col_config:
        df = apply_filter(df, col_name, col_config["filter"])

    return df


def process_boolean(df: pd.DataFrame, col_name: str, col_config: Dict[str, Any]) -> pd.DataFrame:
    """Process boolean column: score, filter."""
    true_value = col_config["true_value"]
    weight = col_config.get("weight", 0.0)

    # Create q-score (weight if true, 0 if false)
    df[f"{col_name}__q"] = (df[col_name] == true_value).astype(float) * weight

    # Apply filter if specified
    if "filter" in col_config:
        df = apply_filter(df, col_name, col_config["filter"])

    return df


def process_size_deviation_source(df: pd.DataFrame, col_name: str, col_config: Dict[str, Any]) -> pd.DataFrame:
    """Process size deviation using pre-calculated cluster_id."""
    raw_deviations = {}

    # Calculate deviations for each cluster
    for cluster_id in df["cluster_id"].unique():
        cluster_genomes = df[df["cluster_id"] == cluster_id]

        if len(cluster_genomes) < 2:
            # Single genome clusters get neutral score
            for idx in cluster_genomes.index:
                raw_deviations[idx] = 0.0
            continue

        # Calculate cluster median size
        cluster_sizes = cluster_genomes[col_name].astype(float)
        median_size = cluster_sizes.median()

        # Calculate deviations
        for idx, row in cluster_genomes.iterrows():
            size_deviation = abs(row[col_name] - median_size) / median_size if median_size > 0 else 0.0
            raw_deviations[idx] = size_deviation

    # Normalize deviations across all clusters
    deviation_values = list(raw_deviations.values())
    if deviation_values and np.std(deviation_values) > 0:
        mean_dev = np.mean(deviation_values)
        std_dev = np.std(deviation_values)

        normalized_deviations = {
            idx: (raw_dev - mean_dev) / std_dev
            for idx, raw_dev in raw_deviations.items()
        }
    else:
        normalized_deviations = {idx: 0.0 for idx in raw_deviations.keys()}

    # Add q-scores only if weight is specified (for scoring)
    weight = col_config.get("weight")
    if weight is not None:
        direction = col_config.get("direction", "lower_better")

        df[f"{col_name}__q"] = pd.Series(normalized_deviations) * weight
        if direction == "higher_better":
            df[f"{col_name}__q"] = -df[f"{col_name}__q"]

        df[f"{col_name}__q"] = df[f"{col_name}__q"].fillna(0.0)

    # Apply filter if specified (filter based on raw deviation, not original column)
    if "filter" in col_config:
        # Create a temporary column with raw deviations for filtering
        df[f"{col_name}__deviation"] = pd.Series(raw_deviations)
        df = apply_filter(df, f"{col_name}__deviation", col_config["filter"])
        # Remove the temporary column
        df = df.drop(columns=[f"{col_name}__deviation"])

    return df


def process_text(df: pd.DataFrame, col_name: str, col_config: Dict[str, Any]) -> pd.DataFrame:
    """Process text column (no scoring, just metadata for plotting)."""
    # Text columns don't contribute to scoring, just pass through
    return df


def process_config(metadata: pd.DataFrame, config: Dict[str, Any],
                   distance_data: pd.DataFrame, target_core_size: int) -> pd.DataFrame:
    """Process configuration step by step with standardized API."""
    df = metadata.copy()

    # Add helper columns
    df["__filter__"] = True
    df["cluster_id"] = calculate_cluster_ids(distance_data, target_core_size).values

    logger.info(f"   Created {target_core_size} phylogenetic clusters")

    # Process each column
    for col_name, col_config in config["columns"].items():
        col_type = col_config.get("type", "text")

        if col_name not in df.columns:
            raise ValueError(f"Column '{col_name}' specified in config but not found in metadata table")

        logger.info(f"   Processing {col_type} column: {col_name}")

        if col_type == "numeric":
            df = process_numeric(df, col_name, col_config)
        elif col_type == "boolean":
            df = process_boolean(df, col_name, col_config)
        elif col_type == "size_deviation_source":
            df = process_size_deviation_source(df, col_name, col_config)
        elif col_type == "text":
            df = process_text(df, col_name, col_config)
        else:
            raise ValueError(f"Unknown column type: {col_type}")

    return df


def calculate_final_scores(df: pd.DataFrame) -> Dict[str, float]:
    """Calculate final quality scores by summing all q-score columns."""
    q_columns = [col for col in df.columns if col.endswith("__q")]

    if not q_columns:
        raise ValueError("No scoring columns found")

    final_scores = df[q_columns].sum(axis=1)
    return final_scores.to_dict()


def select_core_set_api(
    metadata: pd.DataFrame,
    distance_data: pd.DataFrame,
    target_core_size: int,
    config: Dict[str, Any],
) -> Tuple[List[str], Dict[str, float], Dict[str, int], List[str]]:
    """
    Python API to select optimal core set of genomes for OrthoFinder analysis.

    Args:
        metadata: DataFrame with genome metadata
        distance_data: Distance matrix DataFrame with genome IDs as index/columns
        target_core_size: Number of genomes to select for core set
        config: JSON configuration dictionary with columns, weights, filters, plotting

    Returns:
        Tuple of (selected genome IDs, quality scores dict, cluster assignments dict, filtered genome list)
    """
    logger.info(f"🎯 OrthoFinder CoreSelector: Selecting {target_core_size} genomes for core set")

    # Make copies to avoid modifying original data
    metadata = metadata.copy()
    distance_data = distance_data.copy()

    # Validate and normalize configuration
    config = validate_config(config)
    assert metadata.index.name == config['index'], f"{metadata.index.name=} does not match {config['index']=}"

    # Validate sample consistency between metadata and distance matrix
    validate_sample_consistency(metadata, distance_data)

    # Ensure distance matrix matches metadata order
    genome_ids = metadata.index.tolist()
    distance_data = distance_data.loc[genome_ids, genome_ids]

    # Phase 1: Process configuration and calculate scores
    logger.info("\n🏆 Phase 1: Processing configuration and calculating scores...")

    processed_df = process_config(metadata, config, distance_data, target_core_size)

    # Apply filters
    filtered_df = processed_df[processed_df["__filter__"]]
    if len(filtered_df) == 0:
        raise ValueError("No genomes passed the filters")

    filtered_count = len(filtered_df)
    total_count = len(processed_df)
    logger.info(f"   {filtered_count}/{total_count} genomes passed filters")

    # Calculate final quality scores
    quality_scores = calculate_final_scores(filtered_df)

    # Update distance data to match filtered genomes
    filtered_genome_ids = filtered_df.index.tolist()
    filtered_distance_data = distance_data.loc[filtered_genome_ids, filtered_genome_ids]

    # Check if we have enough genomes after filtering
    if len(filtered_df) < target_core_size:
        raise ValueError(f"Only {len(filtered_df)} genomes available after filtering, "
        "but {target_core_size} requested. Adjust filters or reduce target size.")

    # Phase 2: Select genomes using phylogenetic clustering
    logger.info(f"\n🌳 Phase 2: Selecting {target_core_size} genomes using phylogenetic clustering...")

    selected_genomes, cluster_assignments = select_by_phylogenetic_clustering(
        quality_scores=quality_scores,
        distance_data=filtered_distance_data,
        target_core_size=target_core_size,
        metadata=filtered_df
    )

    # Create cluster assignments for ALL genomes (including filtered ones) for visualization
    all_cluster_assignments = {}
    for genome_id in processed_df.index:
        cluster_id = int(processed_df.loc[genome_id, 'cluster_id'])
        all_cluster_assignments[genome_id] = cluster_id

    # Output results
    logger.info(f"\n✅ Selected {len(selected_genomes)} genomes!")

    return selected_genomes, quality_scores, all_cluster_assignments, filtered_df.index.tolist()


def select_core_set(
    distance_matrix: str,
    target_core_size: int,
    output_file: str,
    visualization_file: str,
    metadata_file: Optional[str] = None,
    config_file: Optional[str] = None,
    log_level: str = "INFO",
    max_genomes: Optional[int] = None,
) -> List[str]:
    """
    CLI interface to select optimal core set of genomes for OrthoFinder analysis.

    Supports two modes:
    1. Minimal mode: Only distance matrix - uses random selection within phylogenetic clusters
    2. Quality mode: With metadata and config - uses quality-based selection within clusters

    Args:
        distance_matrix: Path to ANI distance matrix file
        target_core_size: Number of genomes to select for core set
        output_file: File to save selected genome IDs
        visualization_file: Output file for visualization
        metadata_file: Path to TSV file with genome metadata (optional)
        config_file: Path to JSON configuration file (required if metadata_file provided)
        log_level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        max_genomes: Limit analysis to first N genomes (for testing)

    Returns:
        List of selected genome IDs
    """
    # Set up logging
    ortho_logger = logging.getLogger("orthofinder_coreselector")
    ortho_logger.setLevel(log_level.upper())
    if not ortho_logger.handlers:
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(logging.Formatter("%(message)s"))
        ortho_logger.addHandler(handler)

    logger.info(f"📊 Distance matrix: {distance_matrix}")

    # Load distance matrix (always required)
    distance_data = load_distance_matrix(distance_matrix)

    # Check if we have metadata and config for quality-based selection
    if metadata_file and config_file:
        # Quality mode: use metadata and configuration
        logger.info(f"📁 Metadata file: {metadata_file}")
        logger.info(f"⚙️ Config file: {config_file}")

        # Load configuration and metadata
        with open(config_file) as f:
            config = json.load(f)
        config = validate_config(config)
        metadata = load_metadata(metadata_file, index_col=config['index'])

        # Apply max_genomes for debugging (reduce dataset size early)
        if max_genomes:
            logger.debug(f"🔧 Debug mode: Limiting to first {max_genomes} genomes")
            metadata = metadata.head(max_genomes)
            # Get the genome IDs from the reduced metadata
            genome_ids = metadata.index.tolist()
            distance_data = distance_data.loc[genome_ids, genome_ids]

        # Run quality-based selection
        selected_genomes, quality_scores, cluster_assignments, filtered_genome_list = select_core_set_api(
            metadata=metadata,
            distance_data=distance_data,
            target_core_size=target_core_size,
            config=config
        )

    elif metadata_file or config_file:
        # Error: need both metadata and config, or neither
        raise ValueError("If providing metadata_file, config_file is also required. For minimal mode, provide neither.")

    else:
        # Minimal mode: only distance matrix, random selection within clusters
        logger.info("🎲 Minimal mode: Using random selection within phylogenetic clusters")

        # Apply max_genomes for debugging
        if max_genomes:
            logger.debug(f"🔧 Debug mode: Limiting to first {max_genomes} genomes")
            genome_ids = distance_data.index.tolist()[:max_genomes]
            distance_data = distance_data.loc[genome_ids, genome_ids]

        # Create minimal metadata for clustering
        genome_ids = distance_data.index.tolist()
        minimal_metadata = pd.DataFrame(index=genome_ids)

        # Use random quality scores (all equal)
        quality_scores = {genome_id: 1.0 for genome_id in genome_ids}

        # Run phylogenetic clustering with random selection
        selected_genomes, cluster_assignments = select_by_phylogenetic_clustering(
            quality_scores=quality_scores,
            distance_data=distance_data,
            target_core_size=target_core_size,
            metadata=minimal_metadata
        )

        # No filtering in minimal mode
        filtered_genome_list = genome_ids

    # Save selected genomes to output file
    with open(output_file, 'w') as f:
        f.write("\n".join(selected_genomes))

    logger.info(f"\n💾 Results saved to: {output_file}")

    # Create visualization
    if metadata_file and config_file:
        # Full visualization with metadata
        create_visualization(
            metadata=metadata,
            distance_data=distance_data,
            selected_genomes=selected_genomes,
            quality_scores=quality_scores,
            output_file=visualization_file,
            cluster_assignments=cluster_assignments,
            config=config,
            filtered_genomes=filtered_genome_list
        )
    else:
        # Minimal visualization (only phylogenetic tree)
        create_visualization(
            metadata=minimal_metadata,
            distance_data=distance_data,
            selected_genomes=selected_genomes,
            quality_scores=quality_scores,
            output_file=visualization_file,
            cluster_assignments=cluster_assignments,
            config={},  # Empty config for minimal mode
            filtered_genomes=filtered_genome_list
        )

    return selected_genomes


def main():
    fire.Fire(select_core_set)


if __name__ == "__main__":
    main()
