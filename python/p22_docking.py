# -*- coding: utf-8 -*-
"""
p22_docking.py — AutoDock Vina docking of 6 ligands against human Galectin-3 (PDB 1KJL).

Steps (re-runnable, idempotent):
  1. Download ligand SDFs from PubChem PUG-REST (3D preferred, 2D fallback + ETKDG/MMFF).
  2. Prepare receptor: strip HETATM (NAG, GAL, BR, CL, HOH) -> meeko mk_prepare_receptor.
  3. Prepare ligand PDBQTs with meeko mk_prepare_ligand.
  4. Run Vina (box center = NAG heavy-atom centroid from 1KJL, size 25 A, exhaustiveness 16).
  5. Parse poses -> results/62_docking/docking_affinities.csv

Run with: gsmap_env/Scripts/python.exe scripts/p22_docking.py
"""
import csv
import json
import os
import re
import subprocess
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
os.chdir(ROOT)

VINA = ROOT / "tools" / "vina.exe"
SCRIPTS = ROOT / "gsmap_env" / "Scripts"
MK_REC = SCRIPTS / "mk_prepare_receptor.exe"
MK_LIG = SCRIPTS / "mk_prepare_ligand.exe"

PDB_IN = ROOT / "docking" / "1KJL.pdb"
REC_PDB = ROOT / "docking" / "1KJL_protein.pdb"
REC_PDBQT = ROOT / "docking" / "1KJL_protein.pdbqt"

LIG_DIR = ROOT / "docking" / "ligands"
OUT_DIR = ROOT / "results" / "62_docking"
LIG_DIR.mkdir(parents=True, exist_ok=True)
OUT_DIR.mkdir(parents=True, exist_ok=True)

BOX_SIZE = (25, 25, 25)
EXHAUSTIVENESS = 16

# (key, [pubchem names to try in order])
LIGANDS = [
    ("TD139",        ["TD139", "GB0139"]),
    ("lactose",      ["lactose"]),
    ("daunorubicin", ["daunorubicin"]),
    ("dorsomorphin", ["dorsomorphin", "Compound C"]),
    ("mitoxantrone", ["mitoxantrone"]),
    ("triptolide",   ["triptolide"]),
]

PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def urlopen_retry(url, retries=3, timeout=120):
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception as e:
            last = e
            print(f"    [retry {i+1}/{retries}] {url} -> {e}")
    raise last


def compute_box_center():
    xs, ys, zs = [], [], []
    with open(PDB_IN) as f:
        for line in f:
            if line.startswith("HETATM") and line[17:20].strip() == "NAG":
                xs.append(float(line[30:38])); ys.append(float(line[38:46])); zs.append(float(line[46:54]))
    n = len(xs)
    c = (sum(xs)/n, sum(ys)/n, sum(zs)/n)
    print(f"[box] NAG heavy atoms n={n}, center = ({c[0]:.3f}, {c[1]:.3f}, {c[2]:.3f})")
    return c


def prepare_receptor_pdb():
    kept = 0
    with open(PDB_IN) as f, open(REC_PDB, "w") as g:
        for line in f:
            if line.startswith("ATOM") or line.startswith(("TER", "END")):
                g.write(line)
                if line.startswith("ATOM"):
                    kept += 1
    print(f"[receptor] kept {kept} ATOM records -> {REC_PDB}")


def download_ligands():
    """Returns dict key -> {cid, start, sdf_path}"""
    meta = {}
    for key, names in LIGANDS:
        sdf_path = LIG_DIR / f"{key}.sdf"
        got = False
        for name in names:
            # CID lookup
            try:
                cid_data = json.loads(urlopen_retry(f"{PUG}/compound/name/{urllib.parse.quote(name)}/cids/JSON"))
                cid = cid_data["IdentifierList"]["CID"][0]
            except Exception:
                continue
            # 3D first, then 2D
            for rec_type in ("3d", "2d"):
                try:
                    data = urlopen_retry(f"{PUG}/compound/cid/{cid}/SDF?record_type={rec_type}")
                    if data and b"M  END" in data:
                        sdf_path.write_bytes(data)
                        meta[key] = {"cid": cid, "start": rec_type, "name": name, "sdf": str(sdf_path)}
                        print(f"[ligand] {key}: name='{name}' CID={cid} {rec_type} SDF ({len(data)} bytes)")
                        got = True
                        break
                except Exception:
                    continue
            if got:
                break
        if not got:
            print(f"[ligand] {key}: FAILED to download (tried {names})")
    return meta


def ensure_3d(meta):
    """For 2D-start ligands, embed 3D conformer with rdkit ETKDG + MMFF; overwrite SDF."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    for key, m in meta.items():
        if m["start"] == "3d":
            continue
        mol = Chem.SDMolSupplier(m["sdf"], removeHs=False)[0]
        if mol is None:
            mol = Chem.SDMolSupplier(m["sdf"])[0]
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        status = AllChem.EmbedMolecule(mol, params)
        if status != 0:
            raise RuntimeError(f"ETKDG embedding failed for {key}")
        AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
        w = Chem.SDWriter(m["sdf"])
        w.write(mol)
        w.close()
        print(f"[ligand] {key}: 2D->3D embedded (ETKDGv3+MMFF), SDF overwritten")


def prepare_ligands(meta):
    for key, m in meta.items():
        out = LIG_DIR / f"{key}.pdbqt"
        if out.exists() and out.stat().st_size > 0:
            print(f"[meeko] {key}.pdbqt exists, skip")
            continue
        # relative paths: rdkit C++ cannot open non-ASCII (Chinese) absolute paths on Windows
        rel_in = os.path.relpath(m["sdf"], ROOT)
        rel_out = os.path.relpath(out, ROOT)
        r = subprocess.run([str(MK_LIG), "-i", rel_in, "-o", rel_out],
                           capture_output=True, text=True, encoding="utf-8", errors="replace")
        if r.returncode != 0:
            print(f"[meeko] {key} FAILED:\n{r.stdout}\n{r.stderr}")
        else:
            print(f"[meeko] {key}.pdbqt prepared")


def prepare_receptor_pdbqt():
    if REC_PDBQT.exists() and REC_PDBQT.stat().st_size > 0:
        print("[meeko] receptor pdbqt exists, skip")
        return
    r = subprocess.run([str(MK_REC), "--read_pdb", str(REC_PDB), "-o", str(REC_PDBQT.with_suffix("")), "-p",
                        "--default_altloc", "A"],
                       capture_output=True, text=True, encoding="utf-8", errors="replace")
    print(r.stdout)
    if r.returncode != 0:
        print(r.stderr)
        raise RuntimeError("mk_prepare_receptor failed")


def run_vina(meta, center):
    for key in meta:
        ldir = OUT_DIR / key
        ldir.mkdir(exist_ok=True)
        out_pdbqt = ldir / f"{key}_out.pdbqt"
        log = ldir / f"{key}_vina.log"
        cmd = [str(VINA),
               "--receptor", str(REC_PDBQT),
               "--ligand", str(LIG_DIR / f"{key}.pdbqt"),
               "--center_x", f"{center[0]:.3f}", "--center_y", f"{center[1]:.3f}", "--center_z", f"{center[2]:.3f}",
               "--size_x", str(BOX_SIZE[0]), "--size_y", str(BOX_SIZE[1]), "--size_z", str(BOX_SIZE[2]),
               "--exhaustiveness", str(EXHAUSTIVENESS),
               "--out", str(out_pdbqt)]
        r = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace")
        log.write_text(r.stdout + "\n--- STDERR ---\n" + r.stderr, encoding="utf-8")
        if r.returncode != 0:
            print(f"[vina] {key} FAILED (rc={r.returncode}), see {log}")
        else:
            print(f"[vina] {key} done -> {out_pdbqt}")


def parse_results(meta):
    rows = []
    for key, m in meta.items():
        out_pdbqt = OUT_DIR / key / f"{key}_out.pdbqt"
        affs = []
        if out_pdbqt.exists():
            for line in out_pdbqt.read_text().splitlines():
                mm = re.match(r"REMARK VINA RESULT:\s+(-?\d+\.?\d*)", line)
                if mm:
                    affs.append(float(mm.group(1)))
        best = affs[0] if affs else None
        top3 = affs[:3]
        print(f"[result] {key}: best={best}, top3={top3}")
        rows.append({
            "ligand": key,
            "cid": m["cid"],
            "pubchem_name_used": m["name"],
            "structure_start": m["start"],
            "best_affinity_kcal_mol": best,
            "pose2_kcal_mol": top3[1] if len(top3) > 1 else "",
            "pose3_kcal_mol": top3[2] if len(top3) > 2 else "",
            "n_poses": len(affs),
        })
    csv_path = OUT_DIR / "docking_affinities.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"[csv] -> {csv_path}")
    return rows


def main():
    center = compute_box_center()
    prepare_receptor_pdb()
    prepare_receptor_pdbqt()
    meta = download_ligands()
    ensure_3d(meta)
    prepare_ligands(meta)
    run_vina(meta, center)
    rows = parse_results(meta)
    # persist metadata for plotting / report
    with open(OUT_DIR / "docking_meta.json", "w") as f:
        json.dump({"box_center": center, "box_size": BOX_SIZE,
                   "exhaustiveness": EXHAUSTIVENESS, "results": rows}, f, indent=2)


if __name__ == "__main__":
    main()
