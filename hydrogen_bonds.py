from rdkit.Geometry import Point3D
import pandas as pd
import math
from rdkit import Chem
from rdkit.Chem import rdMolTransforms
from collections import defaultdict
import json
from glob import glob
import numpy as np


#load data and save them into json
bonds = defaultdict(list)
for ii, name in enumerate(glob(f"/home/ebart/Dokumenty/preopt-atom-H-soubory/xtb_calculations/*/xtbopt.pdb")):
    protein = Chem.MolFromPDBFile(name, removeHs=False, sanitize=False)
    conf = protein.GetConformer()
    for bond in protein.GetBonds():
        a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()
        if a1.GetSymbol() == "H" or a2.GetSymbol() == "H":
            a1_PDBinfo, a2_PDBinfo = a1.GetPDBResidueInfo(), a2.GetPDBResidueInfo()
            bond_key = "-".join(sorted([f"{a1_PDBinfo.GetResidueName()} {a1_PDBinfo.GetName().strip()}", f"{a2_PDBinfo.GetResidueName()} {a2_PDBinfo.GetName().strip()}"]))
            a1_i,a2_i = a1.GetIdx(), a2.GetIdx()
            bonds[bond_key].append(math.dist(conf.GetAtomPosition(a1_i), conf.GetAtomPosition(a2_i)))
json.dump(bonds, open("bonds_raw.json", 'w'))

with open('bonds_raw.json', 'r') as file:
    data = json.load(file)
    numpy_array = np.array(data)
print(numpy_array)