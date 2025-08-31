"""
Import directly from an Arx folder structure (https://docs.abrinca.com/arx/getting-started/folder-structure-and-metadata/)
"""

import json
import logging
import os
import subprocess
import sys
from pathlib import Path

import fire
import pandas as pd
from arx_tools.folder_looper import FolderGenome, FolderLooper

from orthofinder_coreselector.data_loading import validate_config
from orthofinder_coreselector.main import create_visualization, select_core_set_api

logger = logging.getLogger("orthofinder_coreselector")


DEFAULT_CONFIG = {
    "index": "Identifier",
    "columns": {
        "TaxID": {
            "type": "text",
            "plot": True
        },
        "BUSCO [S]": {
            "weight": 1.0,
            "type": "numeric",
            "normalization": False,
            "direction": "higher_better",
            "filter": {"operator": ">=", "value": 0.85},
            "plot": True,
            "display_unit": "",
            "decimal_places": 3
        },
        "Assembly Nr Scaffolds": {
            "weight": 0.5,
            "type": "numeric",
            "direction": "lower_better",
            "filter": {"operator": "<=", "value": 500},
            "plot": True,
            "decimal_places": 0
        },
        "Representative": {
            "weight": 1.0,
            "type": "boolean",
            "true_value": True,
            "filter": {"operator": "==", "value": True},
            "plot": True
        },
        "Assembly Size": {
            "type": "size_deviation_source",
            "filter": {"operator": "<=", "value": 0.3},
            "plot": True,
            "display_unit": "MB",
            "scale_factor": 1000000,
            "decimal_places": 1
        }
    }
}


def relative_symlink(src, dst):
    src = os.path.abspath(src)
    dst = os.path.abspath(dst)
    relative_path = os.path.relpath(src, os.path.dirname(dst))
    return os.symlink(relative_path, dst)

def absolute_symlink(src, dst):
    src = os.path.abspath(src)
    dst = os.path.abspath(dst)
    return os.symlink(src, dst)


def get_fasta_info(fasta_file:Path) -> (int, int):
    with open(fasta_file) as f:
        fasta = f.read()
    n_scaffolds = fasta.count('>')
    assembly_size = sum([len(line) for line in fasta if not line.startswith('>')])
    return n_scaffolds, assembly_size

def load_metadata(workdir:Path, genome: FolderGenome, link_type: str = 'relative') -> dict:
    if link_type == 'relative':
        link_fn = relative_symlink
    elif link_type == 'absolute':
        link_fn = absolute_symlink
    elif link_type == 'copy':
        import shutil
        link_fn = shutil.copyfile
    else:
        raise ValueError(f"Invalid link type: {link_type}")

    fna = f"{genome.path}/{genome.json['assembly_fasta_file']}"
    faa = f"{genome.path}/{genome.json['cds_tool_faa_file']}"
    link_fn(fna, workdir.joinpath('assemblies').joinpath(genome.identifier))
    link_fn(faa, workdir.joinpath('proteins').joinpath(genome.identifier + '.faa'))
    n_scaffolds, assembly_size = get_fasta_info(fna)
    return {
        'Identifier': genome.identifier,
        'TaxID': genome.organism.json['taxid'],
        'Representative': genome.organism.representative().identifier == genome.identifier,
        'Assembly Size': assembly_size,
        'Assembly Nr Scaffolds': n_scaffolds,
        'BUSCO [S]': genome.json['BUSCO']['S']/genome.json['BUSCO']['T'],
        'Organism': genome.organism.name
    }


def init_orthofinder_arx(
    folder_structure: str | Path,
    outdir: str | Path,
    target_core_size: int,
    max_genomes: int = None,
    link_type: str = 'relative',
    config_file: str | Path | None = None
):
    """
    Import directly from an Arx folder structure (meta_data.tsv and distance_matrix.csv)

    Args:
        folder_structure: Path to the Arx folder structure
        outdir: Output directory for processed data
        target_core_size: Number of genomes to select for core set
        max_genomes: Maximum number of genomes to process (optional)
        link_type: Type of symlinks to create ('relative', 'absolute' or 'copy')
        config_file: Path to custom JSON config file (optional, overrides default config)
    """
    # Convert arguments to Path objects
    folder_structure = Path(folder_structure)
    outdir = Path(outdir)
    if config_file:
        config_file = Path(config_file)

    # Load and validate configuration
    config = DEFAULT_CONFIG.copy()  # Start with default config
    if config_file:
        logger.info(f'📄 Loading custom config from: {config_file}')
        with open(config_file, 'r') as f:
            custom_config = json.load(f)
        # Merge custom config with default (custom overrides default)
        config.update(custom_config)

    # Validate the final configuration
    config = validate_config(config)
    logger.info('✅ Configuration validated successfully')

    logger.setLevel('INFO')
    if not logger.handlers:
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(logging.Formatter("%(message)s"))
        logger.addHandler(handler)

    outdir.mkdir(parents=True, exist_ok=True)
    outdir.joinpath('assemblies').mkdir(exist_ok=True)
    outdir.joinpath('proteins').mkdir(exist_ok=True)
    if not outdir.joinpath('genomes.tsv').exists():
        logger.info('📂 Gathering data and metadata...')
        metadata = []
        folder_looper = FolderLooper(folder_structure)
        for i, genome in enumerate(folder_looper.genomes(skip_ignored=True, sanity_check=False, representatives_only=False), 1):
            print(f"\r{i} {genome.path}\033[K", end='', flush=True)
            metadata.append(load_metadata(outdir, genome, link_type))
            if max_genomes and i >= max_genomes:
                break
        print()  # New line after progress updates
        metadata = pd.DataFrame(metadata)
        metadata.set_index('Identifier', inplace=True)
        metadata.to_csv(outdir.joinpath('genomes.tsv'), sep='\t')
    metadata = pd.read_csv(outdir.joinpath('genomes.tsv'), sep='\t', index_col=0)

    # Generate ANI distance matrix using GenDisCal
    if not outdir.joinpath('ani-distance-matrix.csv').exists():
        logger.info('📊 Creating distance matrix using GenDisCal...')
        subprocess.run(f'GenDisCal --preset approxAni --distancematrix {outdir.joinpath("assemblies")}/* -o {outdir.joinpath("ani-distance-matrix.csv")}', shell=True)
    distance_data = pd.read_csv(outdir.joinpath('ani-distance-matrix.csv'), index_col=0)

    logger.info('🎯 Selecting core set...')
    selected_genomes, quality_scores, cluster_assignments, filtered_genome_list = select_core_set_api(
        metadata = metadata,
        distance_data = distance_data,
        target_core_size = target_core_size,
        config = config,
    )

    # Save selected genomes to output file
    logger.info('💾 Saving selected genomes to output file...')
    with open(outdir.joinpath('selected_genomes.txt'), 'w') as f:
        f.write("\n".join(selected_genomes))

    logger.info("🎨 Visualizing selected genomes on tree...")
    create_visualization(
        metadata, distance_data, selected_genomes, quality_scores, outdir.joinpath('visualization.svg'), cluster_assignments, config, filtered_genome_list
    )



def main():
    fire.Fire(init_orthofinder_arx)


if __name__ == "__main__":
    main()


