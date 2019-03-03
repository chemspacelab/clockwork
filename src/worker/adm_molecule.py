import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem, ChemicalForceFields

import rmsd

def load_sdf(filename):
    """

    load sdf file and return rdkit mol list

    args:
        filename sdf

    return:
        list of mol objs

    """

    suppl = Chem.SDMolSupplier(filename,
                               removeHs=False,
                               sanitize=True)

    mollist = [mol for mol in suppl]

    return mollist


def load_sdf_file(obj):

    suppl = Chem.ForwardSDMolSupplier(obj,
                               removeHs=False,
                               sanitize=True)

    mollist = [mol for mol in suppl]

    return mollist


def save_results(mol, results, filename):

    # make rdkit conformers

    energies = [x[4] for x in results]
    energies = np.array(energies)
    sortidx = np.argsort(energies)

    N_atoms = mol.GetNumAtoms()

    atoms = mol.GetAtoms()
    atoms = [atom.GetSymbol() for atom in atoms]

    coordinates_list = []

    no_hydrogen = np.array([2,5,6])
    no_hydrogen -= 1

    for idx in sortidx:
        result = results[idx]
        coord = result[-1]
        coord = np.array(coord)
        coord = coord.reshape((N_atoms, 3))

        coord -= rmsd.centroid(coord[no_hydrogen])

        coordinates_list.append(coord)

    out = []

    # atoms = np.array(atoms)
    # no_hydrogen = np.where(atoms != 'H')

    for i, coord in enumerate(coordinates_list[1:]):

        U = rmsd.kabsch(coord[no_hydrogen], coordinates_list[0][no_hydrogen])

        coord = np.dot(coord, U)

        strxyz = rmsd.set_coordinates(atoms, coord)
        out += [strxyz]

    f = open(filename, 'w')
    f.write("\n".join(out))
    f.close()

    return

