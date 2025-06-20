# Linux-Bulk-MolecularDocking
Automated Bulk Docking Pipeline
---

This repository provides a automated pipeline for **bulk molecular docking** using [AutoDock Vina](http://vina.scripps.edu/) and optional affinity extraction. It is designed for docking multiple ligands to a single receptor efficiently.

## Project Structure

```bash
.
├── dock.sh                  # ✅ Bulk docking script
├── extract_affinities.sh    # ✅ Extracts best affinities from Vina output
├── gridsize.py              # ✅ Calculates Vina grid box from PDB
├── protein.pdbqt            # 🔁 receptor file in .pdbqt format
├── protein.conf             # 🔁 Vina config file (grid box settings)
├── ligands/                 # 🔁 directory of ligand .pdbqt files
└── docking_results/         # ⚙️ Auto-generated output folder after running `dock.sh`
```

### Legend

* ✅ **Created by this repo** — tools/scripts maintained here
* 🔁 **Provided by user** — inputs that you supply to run the docking
* ⚙️ **Generated** — outputs automatically created during execution


## Setup Instructions

### 1. **Dependencies**

* [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina/releases)
* Python 3.x with `numpy`

### 2. **Prepare Inputs**

* Convert your receptor and ligands to `.pdbqt` format.
* Place all ligand `.pdbqt` files in the `ligands/` folder.
* Create a config file (`protein.conf`) specifying:

  * `center_x`, `center_y`, `center_z`
  * `size_x`, `size_y`, `size_z`

#### Example: `protein.conf`

This file defines the **search space** and **docking settings** for AutoDock Vina. Here's an example:

```conf
center_x = 0.00
center_y = 0.00
center_z = 0.00

size_x = 80.00
size_y = 80.00
size_z = 80.00

energy_range = 3
exhaustiveness = 8
```
Explanation of Fields

| Parameter        | Description                                                                            |
| ---------------- | -------------------------------------------------------------------------------------- |
| `center_x/y/z`   | Coordinates (in Å) for the **center** of the docking box.                              |
| `size_x/y/z`     | Dimensions (in Å) of the docking box along each axis. Larger sizes explore more space. |
| `energy_range`   | Max energy difference (in kcal/mol) between best and worst poses to report.            |
| `exhaustiveness` | Determines **search depth**. Higher values = more accurate but slower.                 |


#### Grid Box Setup

Use the `gridsize.py` script to calculate the optimal center and size values from your receptor `.pdb` file:

```bash
python gridsize.py protein.pdb
```
This prints:

```
CENTER_X = ...
CENTER_Y = ...
...
SIZE_Z = ...
```

Use these in your `protein.conf`.

## Run Bulk Docking

```bash
chmod +x dock.sh
./dock.sh
```

* Each ligand is docked to the receptor using Vina.
* Results are saved in `docking_results/`, with one `.pdbqt` and `.txt` file per ligand.

## Extract Binding Affinities

After docking is complete:

```bash
chmod +x extract_affinities.sh
./extract_affinities.sh
```

* Parses each `*_output.txt` file from Vina
* Extracts the best affinity (lowest energy pose)
* Saves them to `best_affinities.csv` (sorted from strongest to weakest binder)

## Output

You’ll get:

* Docked poses: `ligandname_out.pdbqt`
* Vina stdout: `ligandname_output.txt`
* Affinity summary: `best_affinities.csv`

## Example Output (`best_affinities.csv`)

```csv
Filename,Best Affinity (kcal/mol)
lig1_output.txt,-9.7
lig2_output.txt,-8.5
```
