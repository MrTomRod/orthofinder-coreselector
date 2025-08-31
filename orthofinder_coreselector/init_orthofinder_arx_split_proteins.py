import os

import fire

from .init_orthofinder_arx import relative_symlink


def split_proteins(proteins_folder: str, selected_genomes_file: str, core_folder: str, rest_folder: str):
    # Create target directories
    os.makedirs(core_folder)
    os.makedirs(rest_folder)

    with open(selected_genomes_file) as f:
        selected_genomes = set(f.read().splitlines())

    # loop over all files in proteins_folder
    cnt_core, cnt_rest = 0, 0
    for file in os.listdir(proteins_folder):
        identifier = file.removesuffix('.faa')
        # if file in selected_genomes, link to proteins_core and remove from selected_genomes
        if identifier in selected_genomes:
            relative_symlink(os.path.join(proteins_folder, file), os.path.join(core_folder, file))
            selected_genomes.remove(identifier)
            cnt_core += 1
        else:
            relative_symlink(os.path.join(proteins_folder, file), os.path.join(rest_folder, file))
            cnt_rest += 1

    assert len(selected_genomes) == 0, f"Expected no genomes in selected_genomes, but got {selected_genomes=}"
    print(f"Split {cnt_core} genomes into core and {cnt_rest} genomes into rest")

def main():
    fire.Fire(split_proteins)

if __name__ == "__main__":
    main()
