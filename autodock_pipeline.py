#!/usr/bin/env python3

import os
import sys
import glob
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path

VINA_EXECUTABLE = "./vina"
PROTEIN_FOLDER = "Protein"
LIGAND_FOLDER = "Ligand"
OUTPUT_FOLDER = "docking_results"
CONFIG_FOLDER = "config_files"

EXHAUSTIVENESS = 8
NUM_MODES = 9
ENERGY_RANGE = 3

def parse_pdb(file_path):
    coordinates = []
    with open(file_path, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    coordinates.append((x, y, z))
                except ValueError:
                    continue
    return np.array(coordinates)

def calculate_box_center_and_size(coordinates, padding=5.0):
    if coordinates.size == 0:
        raise ValueError("No valid coordinates found")
    min_coords = np.min(coordinates, axis=0)
    max_coords = np.max(coordinates, axis=0)
    box_size = max_coords - min_coords + padding
    box_center = (max_coords + min_coords) / 2
    return box_center, box_size

def create_config_file(protein_path, output_path, center, size):
    config_content = f"""receptor = {protein_path}
center_x = {center[0]:.3f}
center_y = {center[1]:.3f}
center_z = {center[2]:.3f}
size_x = {size[0]:.3f}
size_y = {size[1]:.3f}
size_z = {size[2]:.3f}
exhaustiveness = {EXHAUSTIVENESS}
num_modes = {NUM_MODES}
energy_range = {ENERGY_RANGE}
"""
    with open(output_path, 'w') as config_file:
        config_file.write(config_content)
    print(f" Config file created: {output_path}")
    return config_content

def extract_best_affinity(log_content):
    try:
        lines = log_content.strip().split('\n')
        affinities = []
        for line in lines:
            line = line.strip()
            if line and line[0].isdigit():
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        affinity = float(parts[1])
                        affinities.append(affinity)
                    except ValueError:
                        continue
        return min(affinities) if affinities else None
    except Exception as e:
        print(f" Error extracting affinity: {e}")
        return None

def run_docking(vina_path, config_path, ligand_path, output_path, log_path):
    cmd = [
        vina_path,
        "--config", config_path,
        "--ligand", ligand_path,
        "--out", output_path
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if result.returncode != 0:
            print(f" WARNING: Docking failed")
            print(f" Error: {result.stderr}")
            return False, None
        log_output = result.stdout + result.stderr
        with open(log_path, 'w') as log_file:
            log_file.write(log_output)
        return True, log_output
    except subprocess.TimeoutExpired:
        print(f" WARNING: Docking timeout")
        return False, None
    except Exception as e:
        print(f" ERROR running docking: {e}")
        return False, None

def main():
    print("="*60)
    print("AUTOMATED MOLECULAR DOCKING PIPELINE")
    print("="*60)
    
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    os.makedirs(CONFIG_FOLDER, exist_ok=True)
    
    if not os.path.exists(VINA_EXECUTABLE):
        print(f"\nERROR: Vina executable not found at: {VINA_EXECUTABLE}")
        print("Please update the VINA_EXECUTABLE path in the script.")
        sys.exit(1)
    
    protein_files = glob.glob(os.path.join(PROTEIN_FOLDER, "*.pdbqt"))
    if not protein_files:
        protein_files = glob.glob(os.path.join(PROTEIN_FOLDER, "*.pdb"))
    if not protein_files:
        print(f"\nERROR: No protein files found in {PROTEIN_FOLDER}")
        print("Please ensure protein files are in PDBQT or PDB format.")
        sys.exit(1)
    
    ligand_files = glob.glob(os.path.join(LIGAND_FOLDER, "*.pdbqt"))
    if not ligand_files:
        print(f"\nERROR: No ligand files found in {LIGAND_FOLDER}")
        print("Please ensure ligand files are in PDBQT format.")
        sys.exit(1)
    
    print(f"\nFound {len(protein_files)} protein(s) and {len(ligand_files)} ligand(s)")
    print(f"Total docking jobs: {len(protein_files) * len(ligand_files)}\n")
    
    all_results = {}
    
    for protein_idx, protein_path in enumerate(protein_files, 1):
        protein_name = Path(protein_path).stem
        print(f"\n[{protein_idx}/{len(protein_files)}] Processing protein: {protein_name}")
        print("-"*60)
        
        protein_output_dir = os.path.join(OUTPUT_FOLDER, protein_name)
        os.makedirs(protein_output_dir, exist_ok=True)
        
        print(" Calculating grid box parameters...")
        try:
            coordinates = parse_pdb(protein_path)
            box_center, box_size = calculate_box_center_and_size(coordinates)
            print(f" Center: ({box_center[0]:.2f}, {box_center[1]:.2f}, {box_center[2]:.2f})")
            print(f" Size: ({box_size[0]:.2f}, {box_size[1]:.2f}, {box_size[2]:.2f})")
        except Exception as e:
            print(f" ERROR: Failed to parse protein file: {e}")
            continue
        
        config_path = os.path.join(CONFIG_FOLDER, f"{protein_name}_config.txt")
        create_config_file(protein_path, config_path, box_center, box_size)
        
        results = []
        successful_docks = 0
        
        for ligand_idx, ligand_path in enumerate(ligand_files, 1):
            ligand_name = Path(ligand_path).stem
            print(f" [{ligand_idx}/{len(ligand_files)}] Docking {ligand_name}...", end=" ")
            
            output_pdbqt = os.path.join(protein_output_dir, f"{ligand_name}_out.pdbqt")
            log_file = os.path.join(protein_output_dir, f"{ligand_name}_log.txt")
            
            success, log_output = run_docking(VINA_EXECUTABLE, config_path, ligand_path,
                                             output_pdbqt, log_file)
            
            if success and log_output:
                best_affinity = extract_best_affinity(log_output)
                if best_affinity is not None:
                    results.append({
                        'Ligand': ligand_name,
                        'Best_Affinity_kcal_mol': best_affinity,
                        'Output_File': output_pdbqt,
                        'Log_File': log_file
                    })
                    successful_docks += 1
                    print(f"✓ (Affinity: {best_affinity:.2f} kcal/mol)")
                else:
                    print("✗ (Failed to extract affinity)")
            else:
                print("✗ (Docking failed)")
        
        if results:
            df = pd.DataFrame(results)
            df = df.sort_values('Best_Affinity_kcal_mol')
            csv_path = os.path.join(protein_output_dir, f"{protein_name}_results.csv")
            df.to_csv(csv_path, index=False)
            all_results[protein_name] = df
            
            print(f"\n Summary for {protein_name}:")
            print(f" - Successful dockings: {successful_docks}/{len(ligand_files)}")
            print(f" - Best ligand: {df.iloc[0]['Ligand']}")
            print(f" - Best affinity: {df.iloc[0]['Best_Affinity_kcal_mol']:.2f} kcal/mol")
            print(f" - Results saved to: {csv_path}")
        else:
            print(f"\n WARNING: No successful dockings for {protein_name}")
    
    print("\n" + "="*60)
    print("FINAL SUMMARY")
    print("="*60)
    
    if all_results:
        summary_data = []
        for protein_name, df in all_results.items():
            best_row = df.iloc[0]
            summary_data.append({
                'Protein': protein_name,
                'Best_Ligand': best_row['Ligand'],
                'Best_Affinity_kcal_mol': best_row['Best_Affinity_kcal_mol'],
                'Total_Ligands_Tested': len(df)
            })
        
        summary_df = pd.DataFrame(summary_data)
        summary_df = summary_df.sort_values('Best_Affinity_kcal_mol')
        summary_path = os.path.join(OUTPUT_FOLDER, "summary_all_proteins.csv")
        summary_df.to_csv(summary_path, index=False)
        
        print(f"\nOverall best results:")
        print(summary_df.to_string(index=False))
        print(f"\nSummary saved to: {summary_path}")
        
        print("\nIndividual protein results available in:")
        for protein_name in all_results.keys():
            print(f" - {OUTPUT_FOLDER}/{protein_name}/{protein_name}_results.csv")
    else:
        print("\nNo successful dockings completed.")
    
    print("\n" + "="*60)
    print("PIPELINE COMPLETED")
    print("="*60)

if __name__ == "__main__":
    main()

