import os
import re
import pandas as pd
import numpy as np
import time

from biopandas.pdb import PandasPdb
from biopandas.mol2 import PandasMol2, split_multimol2

def main():
    start_time = time.time()

    converted_mols = {}
    pdbs = {}
    input_list = pd.read_csv('./data/chemscore2.csv')
    max_distance = 6.0

    mols2 = split_multimol2('./data/chemscore_2_solutions.mol2')

    for pdb in os.listdir('./data/pdbs'):
        name = pdb.split('_')[-1].replace('.pdb', '')
        pdbs[name] = PandasPdb().read_pdb('./data/pdbs/{}'.format(pdb)).df

    for mol2 in mols2:
        pmol = PandasMol2().read_mol2_from_list(mol2_lines=mol2[1], mol2_code=mol2[0])
        converted_mols[pmol.df.iloc[0]['subst_name']] = pmol.df

    def create_atom_dimensional_position(x, y, z):
        return np.array((x, y, z))

    def atom_is_close_to_atom(resiude_atom, atom, distance):
        return np.linalg.norm(resiude_atom - atom) <= distance

    def remove_hydrogen_atoms(df):
        return df['ATOM'][~df['ATOM']['atom_name'].str.contains('H')]

    def get_interactions_molecule_for_residues(molecule, residues):
        matched_atoms = []
        mol_points = []

        for midx, mol_row in molecule.iterrows():
            mol_points.append(create_atom_dimensional_position(
                x=mol_row['x'], y=mol_row['y'], z=mol_row['z'],
            ))

        for pidx, protein_row in remove_hydrogen_atoms(residues).iterrows():
            protein_tag = '{}_{}'.format(
                protein_row['residue_name'], protein_row['residue_number']
            )

            protein_ad = create_atom_dimensional_position(
                x=protein_row['x_coord'], y=protein_row['y_coord'], z=protein_row['z_coord']
            )

            for point in mol_points:
                if atom_is_close_to_atom(protein_ad, point, max_distance) and protein_tag not in matched_atoms:
                    matched_atoms.append(protein_tag)

        return matched_atoms

    def interate_with_expected_residues(item, expected_residues):
        print('> {}'.format(item['NAME']))

        resp = []
        m = converted_mols[item['NAME']]
        r = pdbs[str(item['Gold.Ensemble.ID'])]

        interactions = get_interactions_molecule_for_residues(m, r)

        return [key for key in interactions if key in expected_residues]

    EXPECTED_RESIDUES = ['ALA_98', 'LEU_84', 'ILE_89', 'LEU_91', 'ARG_13', 'ARG_9']
    dataframe = pd.DataFrame(columns=['molecule_name', 'pdb', 'residues', 'residues_quantity'])

    print('Start Molecule Analyze')
    for i, row in input_list.iterrows():
        reactions = interate_with_expected_residues(row, EXPECTED_RESIDUES)
        dataframe = dataframe.append({
            'molecule_name': row['NAME'],
            'pdb': row['Gold.Ensemble.ID'],
            'residues': ", ".join(reactions),
            'residues_quantity': len(reactions),
        }, ignore_index=True)

    body = dataframe.sort_values(by=['residues_quantity'], ascending=False).to_csv(index=False)

    file = open("./data/filter.csv", "w")
    file.write(body)
    file.close()

    print("Process take: %s seconds" % (time.time() - start_time))

if __name__ == '__main__':
    main()