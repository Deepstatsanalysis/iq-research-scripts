import os
import re
import pandas as pd
import numpy as np
import time
import click

from biopandas.pdb import PandasPdb
from biopandas.mol2 import PandasMol2, split_multimol2
from dotenv import load_dotenv

@click.command()
@click.argument('config')
def main(config):
    start_time = time.time()

    env_path = config
    load_dotenv(dotenv_path=env_path)

    converted_mols = {}
    pdbs = {}
    input_list = pd.read_csv(os.getenv('INPUT_LIST'))

    mols2 = split_multimol2(os.getenv('MOL2_FILE'))

    pdb_path = os.getenv('PDBS_FILE_FOLDER')
    for pdb in os.listdir(pdb_path):
        name = pdb.split('_')[-1].replace('.pdb', '')
        pdbs[name] = PandasPdb().read_pdb('{}/{}'.format(pdb_path, pdb)).df

    for mol2 in mols2:
        pmol = PandasMol2().read_mol2_from_list(mol2_lines=mol2[1], mol2_code=mol2[0])
        converted_mols[pmol.df.iloc[0]['subst_name']] = pmol.df

    def create_atom_dimensional_position(x, y, z):
        return np.array((x, y, z))

    def atom_is_close_to_atom(resiude_atom, atom):
        return np.linalg.norm(resiude_atom - atom) <= float(os.getenv('MAX_DISTANCE'))

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
                if atom_is_close_to_atom(protein_ad, point) and protein_tag not in matched_atoms:
                    matched_atoms.append(protein_tag)

        return matched_atoms

    def interate_with_expected_residues(item, expected_residues):
        print('> {}'.format(item['NAME']))

        resp = []
        m = converted_mols[item['NAME']]
        r = pdbs[str(item['Gold.Ensemble.ID'])]

        interactions = get_interactions_molecule_for_residues(m, r)

        return [key for key in interactions if key in expected_residues]

    dataframe = pd.DataFrame(columns=['molecule_name', 'pdb', 'residues', 'residues_quantity'])

    print('Start Molecule Analyze')
    for i, row in input_list.iterrows():
        reactions = interate_with_expected_residues(row, list(os.getenv('RESIDUES').split(",")))
        dataframe = dataframe.append({
            'molecule_name': row['NAME'],
            'pdb': row['Gold.Ensemble.ID'],
            'residues': ", ".join(reactions),
            'residues_quantity': len(reactions),
        }, ignore_index=True)

    body = dataframe.sort_values(by=['residues_quantity'], ascending=False).to_csv(index=False)

    file = open(os.getenv('OUTPUT_LIST'), "w")
    file.write(body)
    file.close()

    print("Process take: %s seconds" % (time.time() - start_time))

if __name__ == '__main__':
    main()