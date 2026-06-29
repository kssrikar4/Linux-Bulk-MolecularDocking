#!/usr/bin/env python3

import os
import glob
import re
import pandas as pd
import numpy as np
from plip.structure.preparation import PDBComplex

DOCKING_RESULTS_FOLDER = "docking_results"
COMPLEXES_FOLDER = "protein_ligand_complexes"
OUTPUT_IMAGES_FOLDER = "interaction_images"
OUTPUT_TABLE = "docking_interactions.xlsx"
ZOOM_BUFFER = 3.0
# Configure visual colors: BG_COLOR can be 'white', 'black', or 'none' (transparent); adjust TEXT_COLOR; PROTEIN_COLOR can be 'default' (multi-color) or a color name.
BG_COLOR = "white"
TEXT_COLOR = "black"
PROTEIN_COLOR = "default"

PDBQT_TO_ELEMENT = {
    'HD': 'H', 'HS': 'H', 'NA': 'N', 'NS': 'N', 
    'OA': 'O', 'OS': 'O', 'SA': 'S'
}

def clean_pdb_for_plip(pdb_path):
    with open(pdb_path, 'r') as f:
        content = f.read()
    cleaned_lines = []
    for line in content.splitlines():
        if line.startswith(('ATOM', 'HETATM')):
            for qttype, element in PDBQT_TO_ELEMENT.items():
                if line.rstrip().endswith(qttype):
                    line = line[:76] + element.rjust(2)
                    break
        cleaned_lines.append(line)
    return '\n'.join(cleaned_lines)

def coord_sel(atom_obj):
    x, y, z = atom_obj.coords
    return f"x > {x-0.02} and x < {x+0.02} and y > {y-0.02} and y < {y+0.02} and z > {z-0.02} and z < {z+0.02}"

def generate_interaction_image(pdb_path, output_path, inter_set):
    import pymol
    from pymol import cmd
    pymol.finish_launching(['pymol', '-qc'])
    cmd.reinitialize()
    cmd.load(pdb_path, "complex")
    
    if BG_COLOR == "none":
        cmd.bg_color("none")
        cmd.set("ray_opaque_background", 0)
    else:
        cmd.bg_color(BG_COLOR)
        cmd.set("ray_opaque_background", 1)
        
    cmd.hide("everything")
    cmd.show("cartoon", "polymer.protein")
    if PROTEIN_COLOR == "default":
        cmd.util.color_chains("polymer.protein")
    elif PROTEIN_COLOR:
        cmd.color(PROTEIN_COLOR, "polymer.protein")
    cmd.set("cartoon_transparency", 0.5)
    cmd.show("sticks", "resn UNL")
    cmd.util.cbag("resn UNL")
    
    hbonds = inter_set.hbonds_ldon + inter_set.hbonds_pdon
    hydrophobics = inter_set.hydrophobic_contacts
    saltbridges = inter_set.saltbridge_lneg + inter_set.saltbridge_pneg
    pistackings = inter_set.pistacking
    pications = inter_set.pication_laro + inter_set.pication_paro
    halogens = inter_set.halogen_bonds
    
    interacting_res = set()
    def add_residue_selection(res_list):
        for r in res_list:
            interacting_res.add((r.restype, r.resnr, r.reschain))
            
    add_residue_selection(hbonds)
    add_residue_selection(hydrophobics)
    add_residue_selection(saltbridges)
    add_residue_selection(pistackings)
    add_residue_selection(pications)
    add_residue_selection(halogens)
    
    selection_parts = []
    for restype, resnr, reschain in interacting_res:
        sel = f"(chain {reschain} and resi {resnr})"
        selection_parts.append(sel)
        cmd.show("sticks", sel)
        cmd.util.cbaw(sel)
        cmd.label(f"{sel} and name CA", f'"{restype}{resnr}"')
        
    cmd.set("label_size", 10)
    cmd.set("label_color", TEXT_COLOR)
    cmd.set("label_font_id", 7)
    
    for idx, hb in enumerate(hbonds):
        dist_name = f"hbond_{idx}"
        cmd.distance(dist_name, coord_sel(hb.a), coord_sel(hb.d))
        cmd.set("dash_color", "orange", dist_name)
        cmd.set("dash_width", 3.0, dist_name)
        cmd.set("dash_gap", 0.15, dist_name)
        
    for idx, hy in enumerate(hydrophobics):
        dist_name = f"hydro_{idx}"
        cmd.distance(dist_name, coord_sel(hy.ligatom), coord_sel(hy.bsatom))
        cmd.set("dash_color", "grey50", dist_name)
        cmd.set("dash_width", 1.5, dist_name)
        cmd.set("dash_gap", 0.2, dist_name)
        cmd.hide("labels", dist_name)
        
    for idx, hal in enumerate(halogens):
        dist_name = f"halogen_{idx}"
        cmd.distance(dist_name, coord_sel(hal.acc), coord_sel(hal.don))
        cmd.set("dash_color", "cyan", dist_name)
        cmd.set("dash_width", 3.0, dist_name)
        cmd.set("dash_gap", 0.15, dist_name)
        
    for idx, sb in enumerate(saltbridges):
        p_name = f"sb_pos_{idx}"
        n_name = f"sb_neg_{idx}"
        dist_name = f"sb_dist_{idx}"
        cmd.pseudoatom(p_name, pos=list(sb.positive.center))
        cmd.pseudoatom(n_name, pos=list(sb.negative.center))
        cmd.distance(dist_name, p_name, n_name)
        cmd.set("dash_color", "red", dist_name)
        cmd.set("dash_width", 3.0, dist_name)
        cmd.hide("everything", p_name)
        cmd.hide("everything", n_name)
        
    for idx, pi in enumerate(pistackings):
        r_name = f"pi_prot_{idx}"
        l_name = f"pi_lig_{idx}"
        dist_name = f"pi_dist_{idx}"
        cmd.pseudoatom(r_name, pos=list(pi.proteinring.center))
        cmd.pseudoatom(l_name, pos=list(pi.ligandring.center))
        cmd.distance(dist_name, r_name, l_name)
        cmd.set("dash_color", "blue", dist_name)
        cmd.set("dash_width", 3.0, dist_name)
        cmd.hide("everything", r_name)
        cmd.hide("everything", l_name)
        cmd.hide("labels", dist_name)
        
    for idx, pic in enumerate(pications):
        r_name = f"pic_ring_{idx}"
        c_name = f"pic_charge_{idx}"
        dist_name = f"pic_dist_{idx}"
        cmd.pseudoatom(r_name, pos=list(pic.ring.center))
        cmd.pseudoatom(c_name, pos=list(pic.charge.center))
        cmd.distance(dist_name, r_name, c_name)
        cmd.set("dash_color", "purple", dist_name)
        cmd.set("dash_width", 3.0, dist_name)
        cmd.hide("everything", r_name)
        cmd.hide("everything", c_name)

    if selection_parts:
        framing_selection = "resn UNL or " + " or ".join(selection_parts)
    else:
        framing_selection = "resn UNL"
        
    cmd.orient("resn UNL")
    cmd.zoom(framing_selection, buffer=ZOOM_BUFFER)
    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_shadows", 1)
    cmd.set("depth_cue", 1)
    cmd.png(output_path, width=1200, height=900, dpi=300, ray=1)

def main():
    os.makedirs(OUTPUT_IMAGES_FOLDER, exist_ok=True)
    complex_files = glob.glob(os.path.join(COMPLEXES_FOLDER, "*_complex.pdb"))
    if not complex_files:
        return
    results = []
    for complex_file in sorted(complex_files):
        filename = os.path.basename(complex_file)
        match = re.match(r"(.+)_(.+)_complex\.pdb", filename)
        if not match:
            continue
        protein_name, ligand_name = match.group(1), match.group(2)
        affinity = None
        with open(complex_file, 'r') as f:
            for line in f:
                if "Binding Affinity:" in line:
                    try:
                        affinity = float(line.split(":")[-1].replace("kcal/mol", "").strip())
                    except ValueError:
                        pass
                    break
        cleaned_pdb_string = clean_pdb_for_plip(complex_file)
        my_complex = PDBComplex()
        try:
            my_complex.load_pdb(cleaned_pdb_string, as_string=True)
            my_complex.analyze()
        except Exception:
            continue
        if not my_complex.interaction_sets:
            results.append({
                'Protein': protein_name, 'Ligand': ligand_name, 'Affinity_kcal_mol': affinity,
                'H-Bonds': 0, 'Hydrophobic': 0, 'Salt_Bridges': 0, 'pi-Stacking': 0,
                'pi-Cation': 0, 'Halogen_Bonds': 0, 'Total_Interacting_Residues': 0,
                'Interacting_Residues': ""
            })
            continue
        lig_key = list(my_complex.interaction_sets.keys())[0]
        inter_set = my_complex.interaction_sets[lig_key]
        num_hbonds = len(inter_set.hbonds_ldon) + len(inter_set.hbonds_pdon)
        num_hydrophobic = len(inter_set.hydrophobic_contacts)
        num_saltbridge = len(inter_set.saltbridge_lneg) + len(inter_set.saltbridge_pneg)
        num_pistacking = len(inter_set.pistacking)
        num_pication = len(inter_set.pication_laro) + len(inter_set.pication_paro)
        num_halogen = len(inter_set.halogen_bonds)
        interacting_res_details = []
        unique_residues = set()
        def add_res(res_list, type_name):
            for r in res_list:
                res_str = f"{r.restype}{r.resnr}"
                unique_residues.add(res_str)
                interacting_res_details.append(f"{res_str}({type_name})")
        add_res(inter_set.hbonds_ldon + inter_set.hbonds_pdon, "H-bond")
        add_res(inter_set.hydrophobic_contacts, "Hydrophobic")
        add_res(inter_set.saltbridge_lneg + inter_set.saltbridge_pneg, "Salt-Bridge")
        add_res(inter_set.pistacking, "pi-Stacking")
        add_res(inter_set.pication_laro + inter_set.pication_paro, "pi-Cation")
        add_res(inter_set.halogen_bonds, "Halogen-Bond")
        interacting_residues_str = ", ".join(sorted(list(set(interacting_res_details))))
        total_interacting = len(unique_residues)
        results.append({
            'Protein': protein_name,
            'Ligand': ligand_name,
            'Affinity_kcal_mol': affinity,
            'H-Bonds': num_hbonds,
            'Hydrophobic': num_hydrophobic,
            'Salt_Bridges': num_saltbridge,
            'pi-Stacking': num_pistacking,
            'pi-Cation': num_pication,
            'Halogen_Bonds': num_halogen,
            'Total_Interacting_Residues': total_interacting,
            'Interacting_Residues': interacting_residues_str
        })
        ligand_dir = os.path.join(OUTPUT_IMAGES_FOLDER, ligand_name)
        os.makedirs(ligand_dir, exist_ok=True)
        image_path = os.path.join(ligand_dir, f"{protein_name}_interaction.png")
        try:
            generate_interaction_image(complex_file, image_path, inter_set)
        except Exception:
            pass
    if results:
        df_out = pd.DataFrame(results)
        with pd.ExcelWriter(OUTPUT_TABLE, engine="openpyxl") as writer:
            for ligand, group in df_out.groupby("Ligand"):
                group = group.sort_values(by=['Protein', 'Affinity_kcal_mol'])
                sheet_name = re.sub(r'[\\/*?:\[\]]', '', ligand)[:31]
                group.to_excel(writer, sheet_name=sheet_name, index=False)

if __name__ == "__main__":
    main()
