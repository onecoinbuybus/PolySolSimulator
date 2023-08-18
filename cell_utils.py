from radonpy.core import utils, poly
import numpy as np
from radonpy.sim import md
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import Geometry as Geom
from copy import deepcopy
import os
import pandas as pd
from scipy.spatial import distance
from scipy.constants import Avogadro
from scipy.spatial.transform import Rotation as R
import random


def get_center_coord_in_cell(mol, cell):
    """
    Adjusts the molecule's coordinates to be at the center of the given cell.

    :param mol: The molecule to be placed at the center.
    :param cell: The cell where the molecule will be placed.
    :return: Adjusted coordinates for the molecule.
    """

    # Compute the center of the cell
    xhi = cell.cell.xhi
    xlo = cell.cell.xlo
    yhi = cell.cell.yhi
    ylo = cell.cell.ylo
    zhi = cell.cell.zhi
    zlo = cell.cell.zlo
    cell_center = [(xhi + xlo) / 2, (yhi + ylo) / 2, (zhi + zlo) / 2]
    # Get the molecule's coordinates and compute its center
    mol_coord = np.array(mol.GetConformer(0).GetPositions())
    mol_center = np.mean(mol_coord, axis=0)

    # Calculate the offset to move the molecule to the cell's center
    offset = cell_center - mol_center
    mol_coord += offset

    return mol_coord


def update_molecule_coordinates(mol, new_coordinates):
    conf = mol.GetConformer(0)
    for i in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(i, new_coordinates[i])


def get_atom_coordinates(mol, confId=0):
    coordinates = []
    conformer = mol.GetConformer(confId)
    for atom in mol.GetAtoms():
        pos = conformer.GetAtomPosition(atom.GetIdx())
        coord = {
            "atom_symbol": atom.GetSymbol(),
            "index": atom.GetIdx(),
            "x": pos.x,
            "y": pos.y,
            "z": pos.z,
        }
        coordinates.append(coord)
    return coordinates


def extract_mol_coordinates(mol, confId=0):
    mol_c = utils.deepcopy_mol(mol)
    mol_c = utils.set_mol_id(mol_c)
    n_mol = utils.count_mols(mol_c)
    coord = np.array(mol_c.GetConformer(confId).GetPositions())

    all_mol_coordinates = {}  # A dictionary to store coordinates for each molecule.

    for i in range(1, n_mol + 1):
        mol_coordinates = []  # List to store coordinates for the current molecule.

        for j, atom in enumerate(mol_c.GetAtoms()):
            if i == atom.GetIntProp("mol_id"):
                mol_coordinates.append(coord[j])

        # Save the coordinates of the current molecule to the dictionary.
        all_mol_coordinates[f"Molecule_{i}"] = mol_coordinates

    return all_mol_coordinates


def find_close_molecules(center_coords, solvent_coords_df, dist_min):
    close_molecules_details = []
    close_molecules = []

    # Ensure center_coords is a NumPy array
    center_coords = np.array(center_coords)

    for col in solvent_coords_df.columns:
        solvent_coords = np.array(solvent_coords_df[col].dropna().tolist())

        for center_atom_index, center_coord in enumerate(center_coords):
            for solvent_atom_index, solvent_atom_coords in enumerate(solvent_coords):
                dist = distance.euclidean(center_coord, solvent_atom_coords)
                if dist < dist_min:
                    details = {
                        "molecule": col,
                        "center_atom_index": center_atom_index,
                        "solvent_atom_index": solvent_atom_index,
                        "distance": dist,
                    }
                    close_molecules_details.append(details)
                    close_molecules.append(col)
                    break

    # Removing duplicates from close_molecules
    close_molecules = list(set(close_molecules))
    return close_molecules, close_molecules_details


def remove_close_molecules(mol, close_molecules):
    mol_c = Chem.RWMol(mol)
    mol_c = utils.set_mol_id(mol_c)
    atoms_to_remove = []

    for molecule in close_molecules:
        mol_id = int(molecule.split("_")[1])

        for atom in mol_c.GetAtoms():
            if atom.GetIntProp("mol_id") == mol_id:
                atoms_to_remove.append(atom.GetIdx())

    for idx in sorted(atoms_to_remove, reverse=True):
        mol_c.RemoveAtom(idx)

    return mol_c.GetMol()


def randomly_delete_molecules(mol, num_molecules_to_delete):
    max_mol_id = max([atom.GetIntProp("mol_id") for atom in mol.GetAtoms()])
    ids_to_delete = random.sample(range(1, max_mol_id + 1), num_molecules_to_delete)
    atom_indices_to_delete = []
    for atom in mol.GetAtoms():
        if atom.GetIntProp("mol_id") in ids_to_delete:
            atom_indices_to_delete.append(atom.GetIdx())

    rw_mol = Chem.RWMol(mol)

    for idx in sorted(atom_indices_to_delete, reverse=True):
        rw_mol.RemoveAtom(idx)
    return rw_mol


def combine_mols(mol1, mol2, res_name_1="RU0", res_name_2="RU0"):
    """
    poly.combine_mols

    Combining mol1 and mol2 taking over the angles, dihedrals, impropers, and cell data

    Args:
        mol1, mol2: RDkit Mol object

    Optional args:
        res_name_1, res_name_2: Set residue name of PDB

    Returns:
        RDkit Mol object
    """

    mol = Chem.rdmolops.CombineMols(mol1, mol2)

    mol1_n = mol1.GetNumAtoms()
    mol2_n = mol2.GetNumAtoms()
    angles = []
    dihedrals = []
    impropers = []
    cell = None

    if hasattr(mol1, "angles"):
        angles = deepcopy(mol1.angles)
    if hasattr(mol2, "angles"):
        for angle in mol2.angles:
            angles.append(
                utils.Angle(
                    a=angle.a + mol1_n,
                    b=angle.b + mol1_n,
                    c=angle.c + mol1_n,
                    ff=deepcopy(angle.ff),
                )
            )

    if hasattr(mol1, "dihedrals"):
        dihedrals = deepcopy(mol1.dihedrals)
    if hasattr(mol2, "dihedrals"):
        for dihedral in mol2.dihedrals:
            dihedrals.append(
                utils.Dihedral(
                    a=dihedral.a + mol1_n,
                    b=dihedral.b + mol1_n,
                    c=dihedral.c + mol1_n,
                    d=dihedral.d + mol1_n,
                    ff=deepcopy(dihedral.ff),
                )
            )

    if hasattr(mol1, "impropers"):
        impropers = deepcopy(mol1.impropers)
    if hasattr(mol2, "impropers"):
        for improper in mol2.impropers:
            impropers.append(
                utils.Improper(
                    a=improper.a + mol1_n,
                    b=improper.b + mol1_n,
                    c=improper.c + mol1_n,
                    d=improper.d + mol1_n,
                    ff=deepcopy(improper.ff),
                )
            )

    if hasattr(mol1, "cell"):
        cell = deepcopy(mol1.cell)
    elif hasattr(mol2, "cell"):
        cell = deepcopy(mol2.cell)

    # Generate PDB information and repeating unit information
    resid = []
    for i in range(mol1_n):
        atom = mol.GetAtomWithIdx(i)
        atom_name = (
            atom.GetProp("ff_type") if atom.HasProp("ff_type") else atom.GetSymbol()
        )
        if atom.GetPDBResidueInfo() is None:
            atom.SetMonomerInfo(
                Chem.AtomPDBResidueInfo(
                    atom_name,
                    residueName=res_name_1,
                    residueNumber=1,
                    isHeteroAtom=False,
                )
            )
            resid.append(1)
        else:
            atom.GetPDBResidueInfo().SetName(atom_name)
            resid1 = atom.GetPDBResidueInfo().GetResidueNumber()
            resid.append(resid1)

    max_resid = max(resid) if len(resid) > 0 else 0
    for i in range(mol2_n):
        atom = mol.GetAtomWithIdx(i + mol1_n)
        atom_name = (
            atom.GetProp("ff_type") if atom.HasProp("ff_type") else atom.GetSymbol()
        )
        if atom.GetPDBResidueInfo() is None:
            atom.SetMonomerInfo(
                Chem.AtomPDBResidueInfo(
                    atom_name,
                    residueName=res_name_2,
                    residueNumber=1 + max_resid,
                    isHeteroAtom=False,
                )
            )
            resid.append(1 + max_resid)
        else:
            atom.GetPDBResidueInfo().SetName(atom_name)
            resid2 = atom.GetPDBResidueInfo().GetResidueNumber()
            atom.GetPDBResidueInfo().SetResidueNumber(resid2 + max_resid)
            resid.append(resid2 + max_resid)

    max_resid = max(resid) if len(resid) > 0 else 0
    mol.SetIntProp("num_units", max_resid)

    setattr(mol, "angles", angles)
    setattr(mol, "dihedrals", dihedrals)
    setattr(mol, "impropers", impropers)
    if cell is not None:
        setattr(mol, "cell", cell)

    return mol


def simple_cell(mol, n, density=0.03):
    """
    Simple unit cell generator for single component system based on density

    Args:
        mol: RDkit Mol object
        n: Number of molecules in the unit cell (int)
        density: (float, g/cm3)

    Returns:
        Rdkit Mol object
    """

    # Calculate the mass of the molecule
    mol_mass = sum([atom.GetMass() for atom in mol.GetAtoms()])
    mol_volume = mol_mass / (density * Avogadro) * 1e24  # in Angstrom^3
    edge_length = mol_volume ** (1 / 3)
    grid_size = int(n ** (1 / 3))

    cell_c = Chem.Mol()
    xhi, yhi, zhi = (
        edge_length * grid_size,
        edge_length * grid_size,
        edge_length * grid_size,
    )
    xlo, ylo, zlo = 0, 0, 0

    # Assuming utils.Cell is a project-specific class
    setattr(cell_c, "cell", utils.Cell(xhi, xlo, yhi, ylo, zhi, zlo))

    for i in range(grid_size):
        for j in range(grid_size):
            for k in range(grid_size):
                mol_coord = np.array(mol.GetConformer(0).GetPositions())

                # Random rotation
                rot = R.from_euler("xyz", np.random.uniform(-np.pi, np.pi, 3))
                mol_coord = rot.apply(mol_coord)

                # Translation to the grid position
                trans = np.array([i * edge_length, j * edge_length, k * edge_length])
                mol_coord += trans

                cell_n = cell_c.GetNumAtoms()

                # Assuming combine_mols is a project-specific function
                cell_c = combine_mols(cell_c, mol)

                # Set atomic coordinate
                for l in range(mol.GetNumAtoms()):
                    cell_c.GetConformer(0).SetAtomPosition(
                        cell_n + l,
                        Geom.Point3D(mol_coord[l, 0], mol_coord[l, 1], mol_coord[l, 2]),
                    )

    return cell_c


def create_mixture_cell(polymer_mol, solvent_cell, min_dist_threshold=3.0):
    """
    This function takes a polymer molecule and solvent cell as input,
    and returns a mixture cell. It handles moving the polymer, extracting
    solvent coordinates, identifying close molecules based on a minimum distance
    threshold, removing those close molecules, and finally combining the moved
    polymer with the updated cell.

    Parameters:
    polymer_mol: The polymer molecule.
    solvent_cell: The solvent cell.
    min_dist_threshold: The minimum distance threshold for identifying close molecules.

    Returns:
    mixture_cell: The combined cell of polymer and solvent.
    """
    cell_copy = utils.deepcopy_mol(solvent_cell)
    center_coords = get_center_coord_in_cell(polymer_mol, cell_copy)
    moved_polymer = utils.deepcopy(polymer_mol)
    update_molecule_coordinates(moved_polymer, center_coords)
    solvent_coords = extract_mol_coordinates(cell_copy)
    solvent_coords_df = pd.DataFrame(solvent_coords)
    close_mols, _ = find_close_molecules(
        center_coords, solvent_coords_df, min_dist_threshold
    )
    updated_cell = remove_close_molecules(cell_copy, close_mols)
    mixture_cell = combine_mols(moved_polymer, updated_cell)

    return mixture_cell

    """
    # Example:
    # from radonpy.core import utils
    # from rdkit import Chem
    # import cell_utils
    # solvent = 'CC(=O)C'
    # mol_solvent = utils.mol_from_smiles(solvent)
    # homopoly = Chem.MolFromPDBFile('homopoly.pdb')
    # solvent_cell = cell_utils.simple_cell(mol_solvent, 500, density=0.2)
    # mixture_cell = create_mixture_cell(polymer_mol=homopoly, solvent_cell=solvent_cell, min_dist_threshold=10.0)
    """
