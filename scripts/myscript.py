import os
import pandas as pd

mol_csv = pd.read_csv("C:/Users/Pilsen/Desktop/projeto_aedes/estruturas_para_sintese/resultados/corrigido/solutions_chemscore2_corr.csv", names=['Column5']).Column5.tolist()
mol_list = [name.replace(".mol2", "") for name in os.listdir("C:/Users/Pilsen/Desktop/projeto_aedes/estruturas_para_sintese/3D/mopac/mopac_out")]
not_found = []

# for mol in mol_list:
#     if mol not in mol_csv:
#         not_found.append(mol)

for mol in mol_csv:
    if mol not in mol_list:
        not_found.append(mol)

print(not_found)