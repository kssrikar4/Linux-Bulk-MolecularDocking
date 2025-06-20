#!/bin/bash

receptor="protein.pdbqt"
config="protein.conf"
output_dir="docking_results"

if ! ls ligands/*.pdbqt 1> /dev/null 2>&1; then
    echo "No ligand .pdbqt files found in ligands/ directory. Exiting."
    exit 1
fi

mkdir -p "$output_dir"

for ligand in ligands/*.pdbqt; do
    ligand_name=$(basename "$ligand" .pdbqt)
    output_file="${output_dir}/${ligand_name}_out.pdbqt"
    log_file="${output_dir}/${ligand_name}_output.txt"

    vina_path="path/to/vina_1.2.7_linux_x86_64"  #Change with the path of vina in your system
    $vina_path --receptor "$receptor" --ligand "$ligand" --config "$config" --out "$output_file" > "$log_file"

    if [[ $? -ne 0 ]]; then
        echo "Docking failed for $ligand_name. Skipping."
        continue
    fi

    echo "Docking completed for $ligand_name"
done

echo "Mass docking completed. Results are in $output_dir."

