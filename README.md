# Automated Bulk Molecular Docking Pipeline

A automated pipeline for bulk molecular docking using AutoDock Vina. This toolkit enables high-throughput docking of multiple ligands against multiple receptor targets with minimal user intervention.

## Project Structure

```
.
├── Protein/                    # Receptor files (.pdbqt or .pdb)
├── Ligand/                     # Ligand files (.pdbqt)
├── config_files/               # Auto-generated configuration files
├── docking_results/            # Docking outputs and analysis
├── protein_ligand_complexes/   # Best protein-ligand complexes
├── interaction_images/         # PyMOL interaction images
├── gridsize.py                 # Grid parameter calculator
├── autodock_pipeline.py        # Main docking pipeline
├── combine_complexes.py        # Complex structure builder
├── analyze_interactions.py     # Interaction analyzer and visualizer
├── extract_affinities.sh       # Affinity data extractor
└── docking_interactions.xlsx    # Interaction profiling Excel report
```

## Setup Instructions

### 1. Prerequisites

- [**AutoDock Vina**](https://github.com/ccsb-scripps/AutoDock-Vina/releases)
- **Python 3.7+** with packages: `numpy`, `pandas`, `openpyxl`, `plip`
- [**PyMOL**](https://pymol.org/) (for interaction visualization)
- **Open Babel** (for file format conversion)

### 2. Input Preparation

#### Convert Ligands to PDBQT Format
```bash
# Install Open Babel
sudo apt update
sudo apt install -y openbabel python3-openbabel pymol
```

- Convert all SDF files in Ligand folder
```bash
for f in *.sdf; do 
    echo "Converting $f..."; 
    obabel "$f" -O "${f%.sdf}.pdbqt"; 
done
```

or 

- Convert all PDB files in Ligand folder
```bash
for f in *.pdb; do 
    echo "Converting $f..."; 
    obabel "$f" -O "${f%.pdb}.pdbqt"; 
done
```

#### Prepare Receptor Files
```bash
# Install AutoDockTools (one-time setup)
python -m pip install git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3

# Convert receptor PDB files to PDBQT
for f in *.pdb; do
    echo "Preparing ${f%.pdb}.pdbqt ..."
    python -m AutoDockTools.Utilities24.prepare_receptor4 -r "$f" -o "${f%.pdb}.pdbqt" -A checkhydrogens
done
```

### 3. Run Docking Pipeline

Execute the main automated docking workflow:
```bash
#one-time setup
pip install numpy pandas openpyxl plip pymol-open-source

python autodock_pipeline.py
```

### 4. Generate Complex Structures

After docking completion, create protein-ligand complex files:
```bash
python combine_complexes.py
```

### 5. Analyze Interactions and Visualizations

Analyze interaction profiles and generate PyMOL-rendered 3D visualization images:
```bash
python analyze_interactions.py
```

## Output Description

### Primary Results
- **docking_results/**: Individual docking outputs organized by protein
  - PDBQT files of docked poses
  - Detailed docking logs
  - CSV files with affinity rankings
  - Summary CSV across all proteins

### Complex Structures  
- **protein_ligand_complexes/**: Ready-to-visualize structures
  - Best-scoring protein-ligand complexes in PDB format
  - Suitable for molecular visualization software

### Interaction Analysis
- **interaction_images/**: PyMOL-rendered 3D visualization of protein-ligand interactions organized by ligand
- **docking_interactions.xlsx**: Multi-sheet Excel workbook compiling interaction profiles, bond/contact counts, and list of interacting residues for each ligand

## Utility Scripts

- **gridsize.py**: Calculate optimal grid box parameters for docking
- **extract_affinities.sh**: Extract and sort binding affinities from log files

## Features

- **High-Throughput**: Process multiple proteins and ligands simultaneously
- **Automated Grid Detection**: Automatic calculation of binding site boxes
- **Comprehensive Reporting**: Ranked results with detailed affinity data
- **Complex Generation**: Ready-to-use structures for visualization
- **Interaction Profiling**: Automated detection of hydrogen bonds, hydrophobic contacts, salt bridges, pi-interactions, and halogen bonds, paired with PyMOL-rendered 3D interaction diagrams
- **Error Resilience**: Robust error handling and progress tracking

## Usage Notes

- Ensure Vina executable is in your PATH or update `VINA_EXECUTABLE` in the script
- Protein and ligand files should be properly prepared (charges, hydrogens)
- Default docking parameters can be modified in `autodock_pipeline.py`

## Acknowledgements

This pipeline utilizes:

- **AutoDock Vina** - For molecular docking calculations ([DOI: 10.1002/jcc.21334](https://doi.org/10.1002/jcc.21334))
- **AutoDockTools_py3** - For PDB to PDBQT file preparation ([GitHub](https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3))
- **PLIP (Protein-Ligand Interaction Profiler)** - For interaction analysis ([GitHub](https://github.com/pharmai/plip))
- **PyMOL** - For molecular visualization and figure generation ([Website](https://pymol.org/))

## License

This project is provided for academic and research use. Please cite the original tools (Vina and AutoDockTools) in any publications resulting from their use.
