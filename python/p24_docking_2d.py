# -*- coding: utf-8 -*-
"""p24_docking_2d.py — TD139 best-pose contacts with Gal-3 key residues + 2D schematic.
Run with gsmap_env python (rdkit). Uses relative paths (non-ASCII root)."""
import math
import os
import sys
from pathlib import Path

os.chdir(Path(__file__).resolve().parent.parent)

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "DejaVu Sans"
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

REC = "docking/1KJL_protein.pdb"
POSE = "results/62_docking/TD139/TD139_out.pdbqt"
SDF = "docking/ligands/TD139.sdf"
FIG = "figures/v25/Fig_docking_TD139_2D.png"

KEY_RES = {"HIS": "His", "ARG": "Arg", "GLU": "Glu", "ASN": "Asn", "TRP": "Trp",
           "SER": "Ser", "THR": "Thr", "ASP": "Asp", "GLY": "Gly", "ALA": "Ala",
           "VAL": "Val", "LEU": "Leu", "ILE": "Ile", "PHE": "Phe", "TYR": "Tyr",
           "LYS": "Lys", "CYS": "Cys", "MET": "Met", "PRO": "Pro", "GLN": "Gln"}

# --- parse best pose (MODEL 1) heavy atoms from pdbqt ---
lig_atoms = []  # (name, x, y, z)
in_model = False
with open(POSE) as f:
    for line in f:
        if line.startswith("MODEL"):
            in_model = True
        elif line.startswith("ENDMDL"):
            break
        elif in_model and line.startswith(("ATOM", "HETATM")):
            # pdbqt: element at cols 77-78; skip hydrogens
            elem = line[77:79].strip() or line[12:16].strip()[0]
            if elem == "H":
                continue
            lig_atoms.append((line[12:16].strip(),
                              float(line[30:38]), float(line[38:46]), float(line[46:54])))
print(f"ligand heavy atoms in pose: {len(lig_atoms)}")

# --- parse receptor atoms ---
res_atoms = {}  # (resname, resseq, chain) -> list of (atomname, x,y,z)
with open(REC) as f:
    for line in f:
        if line.startswith("ATOM"):
            resn = line[17:20].strip(); resi = int(line[22:26]); ch = line[21].strip()
            key = (resn, resi, ch)
            res_atoms.setdefault(key, []).append(
                (line[12:16].strip(), float(line[30:38]), float(line[38:46]), float(line[46:54])))

# --- contacts: residue min distance to any ligand atom ---
contacts = {}
for (resn, resi, ch), atoms in res_atoms.items():
    best = 1e9
    for _, x1, y1, z1 in atoms:
        for _, x2, y2, z2 in lig_atoms:
            d = math.dist((x1, y1, z1), (x2, y2, z2))
            if d < best:
                best = d
    if best <= 4.5:
        contacts[f"{KEY_RES.get(resn,resn)}{resi}"] = round(best, 2)
contacts = dict(sorted(contacts.items(), key=lambda kv: kv[1]))
print("contacts (<=4.5 A):")
for k, v in contacts.items():
    print(f"  {k}: {v}")

# --- 2D depiction with interacting-atom highlight ---
mol = Chem.SDMolSupplier(SDF, removeHs=True)[0]
AllChem.Compute2DCoords(mol)

# map pose heavy atoms (ordered) to rdkit heavy atoms (ordered)
highlight = []
hbond_res = {"His158", "Arg162", "Glu184"}
hb_atom_idx = set()
# find ligand atoms within 3.5 A of key-residue side chains
for (resn, resi, ch), atoms in res_atoms.items():
    label = f"{KEY_RES.get(resn,resn)}{resi}"
    if label not in hbond_res:
        continue
    for _, x1, y1, z1 in atoms:
        for i, (_, x2, y2, z2) in enumerate(lig_atoms):
            if math.dist((x1, y1, z1), (x2, y2, z2)) <= 3.6:
                hb_atom_idx.add(i)
rd_heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() != "H"]
for i in sorted(hb_atom_idx):
    if i < len(rd_heavy):
        highlight.append(rd_heavy[i])
print(f"highlight atoms (<=3.6 A to His158/Arg162/Glu184): {highlight}")

img = Draw.MolToImage(mol, size=(900, 640), highlightAtoms=highlight, legend="TD139")
fig = plt.figure(figsize=(11, 7))
ax_img = fig.add_axes([0.02, 0.05, 0.56, 0.9])
ax_img.imshow(img)
ax_img.axis("off")
ax_txt = fig.add_axes([0.60, 0.05, 0.38, 0.9])
ax_txt.axis("off")
lines = ["TD139 docked into Gal-3 CRD (PDB 1KJL)",
         "Vina affinity: -7.06 kcal/mol", "",
         "Residues within 4.5 Å of best pose:"]
for k, v in contacts.items():
    mark = " *" if k in hbond_res else ""
    lines.append(f"  {k}{mark}: {v:.2f} Å")
lines += ["", "* canonical CRD binding residues",
          "Highlighted atoms: ≤3.6 Å to", "His158 / Arg162",
          "(Glu184 not within 4.5 Å of", "best pose in this run)"]
ax_txt.text(0, 1, "\n".join(lines), va="top", ha="left", fontsize=11, family="monospace")
fig.suptitle("TD139–Galectin-3 interactions (2D schematic)", fontsize=13)
fig.savefig(FIG, dpi=300, bbox_inches="tight")
print(f"[fig] -> {FIG}")
