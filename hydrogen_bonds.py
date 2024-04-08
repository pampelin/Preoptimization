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

###Zde vypisuju co mám v JSON
#with open('bonds_raw.json', 'r') as file:
#    data = json.load(file)
#    numpy_array = np.array(data)
#print(numpy_array)
###Konec vypisování co mám v JSON


with open('bonds_raw.json', 'r') as file:
    bonds_data = json.load(file)

# Výpočet průměrných délek vazeb
average_bonds_lengths = {}
for bond_key, lengths in bonds_data.items():
    average_length = np.mean(lengths)
    average_bonds_lengths[bond_key] = average_length

# Uložení průměrů do nového JSON souboru
output_path = 'average_bonds_lengths.json'
with open(output_path, 'w') as file:
    json.dump(average_bonds_lengths, file)

print(f'Průměrné délky vazeb byly uloženy do souboru: {output_path}')


###Grafy, které se mi nedaří vykreslit
#file_path = 'bonds_raw.json'
#df = pd.read_json(file_path)

#print(df.head())

#x_column = 'x_column_name'
#y_column = 'y_column_name'

#plt.figure(figsize=(10, 6))
#plt.plot(df[x_column], df[y_column], marker='o')
#plt.title('Graf - vazby')
#plt.xlabel(x_column)
#plt.ylabel(y_column)
#plt.show()

##další zkouška
#with open('bonds_raw.json', 'r') as file:
#    data = json.load(file)
#    numpy_array = np.array(data)
#print(numpy_array)

#dictionary = json.load(open('bonds_raw.json', 'r'))
#xAxis = [key for key, value in dictionary.items()]
#yAxis = [value for key, value in dictionary.items()]
#plt.grid(True)
#
### LINE GRAPH ##
#plt.plot(xAxis,yAxis, color='maroon', marker='o')
#plt.xlabel('variable')
#plt.ylabel('value')
#
#plt.show()

