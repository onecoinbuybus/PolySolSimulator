from radonpy.core import utils, poly, calc
from radonpy.ff.gaff2_mod import GAFF2_mod
from radonpy.sim import qm
import numpy as np
from collections import Counter
from collections import defaultdict
from rdkit import Chem
import cell_utils

ff = GAFF2_mod()

smiles = "*C(C*)c1ccccc1"
ter_smiles = "*C"

temp = 300
press = 1.0
omp_psi4 = 20
mpi = 20
omp = 20
gpu = 0
mem = 100000
chain_length = 5

from polysol_md import Poly_Sol_EMD, make_work_dir
from charge_utils import AtomIndexing, set_charges, get_selected_charges_list

main_work_path = make_work_dir("test_run")
print("main_work_path", main_work_path)

work_dir = main_work_path


mol = utils.mol_from_smiles(smiles)

mol, energy = qm.conformation_search(
    mol,
    ff=ff,
    work_dir=work_dir,
    psi4_omp=omp_psi4,
    mpi=mpi,
    omp=omp,
    memory=mem,
    log_name="mol",
    etkdg_ver=2,
)

mol1 = utils.deepcopy(mol)

atom_indexing = AtomIndexing(mol)
sorted_result_dict, processed_mol, homopoly_rw = atom_indexing.process_mol()

qm.assign_charges(
    homopoly_rw,
    charge="gasteiger",
    opt=False,
    work_dir=work_dir,
    omp=omp_psi4,
    memory=mem,
    log_name="monomer1",
)


selected_indexes = list(sorted_result_dict.values())
selected_charges = get_selected_charges_list(
    homopoly_rw, selected_indexes, charge="gasteiger"
)

mol1 = set_charges(mol1, selected_charges, charge_type="gasteiger")
utils.MolToJSON(mol1, "ps.json")
print("finish")


solvent = "CC(=O)C"
mol_solvent = utils.mol_from_smiles(solvent)
mol_solvent, energy = qm.conformation_search(
    mol_solvent,
    ff=ff,
    work_dir=work_dir,
    psi4_omp=omp_psi4,
    mpi=mpi,
    omp=omp,
    memory=mem,
    log_name="mol_solvent",
    etkdg_ver=2,
)
qm.assign_charges(
    mol_solvent,
    charge="gasteiger",
    opt=False,
    work_dir=work_dir,
    omp=omp_psi4,
    memory=mem,
    log_name="solvent1",
)


result_solvent = ff.ff_assign(mol_solvent)
# solvent_cell = cell_utils.simple_cell(mol_solvent, 1000, density=0.05)

solvent_cell = poly.amorphous_cell(mol_solvent, 1000, density=0.05)


# Electronic propety calculation  Mulliken
qm.assign_charges(
    mol1,
    charge="gasteiger",
    opt=False,
    work_dir=work_dir,
    omp=omp_psi4,
    memory=mem,
    log_name="monomer1",
)
# RESP charge calculation of a termination unit
ter = utils.mol_from_smiles(ter_smiles)
qm.assign_charges(
    ter,
    charge="gasteiger",
    opt=True,
    work_dir=work_dir,
    omp=omp_psi4,
    memory=mem,
    log_name="ter1",
)
polymer_mol = poly.polymerize_rw(mol1, 20, tacticity="atactic")
polymer_mol = poly.terminate_rw(polymer_mol, ter)
result = ff.ff_assign(polymer_mol)


poly_sol_emd = Poly_Sol_EMD(
    solvent_cell, polymer_mol, main_work_path, ff, omp_psi4, mpi, omp, mem, gpu
)


protocol = "md1"
min_dist_threshold = 5
result = poly_sol_emd.exec(protocol, min_dist_threshold)

print("finish all emd")
