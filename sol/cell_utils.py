import numpy as np
from scipy.constants import Avogadro
from rdkit import Chem, Geometry as Geom
from scipy.spatial.transform import Rotation as R
from radonpy import utils


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
