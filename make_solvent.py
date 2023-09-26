import os
import yaml
import argparse
from radonpy.core import utils, poly
from radonpy.ff.gaff2_mod import GAFF2_mod
from radonpy.sim import qm
from polysol_md import SolventEquilibration, make_work_dir
import cell_utils
from datetime import datetime

def main(config_file, solvent_smiles):
    with open(config_file, "r") as file:
        config = yaml.safe_load(file)
        
    ff = GAFF2_mod()
    omp_psi4 = config["omp_psi4"]
    mpi = config["mpi"]
    omp = config["omp"]
    gpu = config["gpu"]
    mem = config["mem"]
    num_molecules = config["num_molecules"]
    density = config["density"]
    amorphous_cell = config["amorphous_cell"]
    simple_cell = config["simple_cell"]

    main_work_path = make_work_dir(config["main_work_dir_name"])
    sol_md_path = os.path.join(main_work_path, "sol_md")
    work_dir = sol_md_path
    start_time = datetime.now()
    print('start make solvent', start_time)

    mol_solvent = utils.mol_from_smiles(solvent_smiles)
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
    #gasteiger
    qm.assign_charges(
        mol_solvent,
        charge="RESP",
        opt=False,
        work_dir=work_dir,
        omp=omp_psi4,
        memory=mem,
        log_name="mol_solvent",
    )
    result_solvent = ff.ff_assign(mol_solvent)

    if amorphous_cell:
        solvent_cell = poly.amorphous_cell(mol_solvent, num_molecules, density=density)
    elif simple_cell:
        solvent_cell = cell_utils.simple_cell(
            mol_solvent, num_molecules, density=density
        )

    sol_equil = SolventEquilibration(
        solvent_cell, None, work_dir, ff, omp_psi4, mpi, omp, mem, gpu
    )
    result = sol_equil.exec()
    utils.MolToJSON(result, "solvent.json")

    print("finish make solvent, time cost:", datetime.now() - start_time)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Solvent Equilibration")
    parser.add_argument("config_file", help="Path to the YAML configuration file")
    parser.add_argument("solvent_smiles", help="SMILES string of the solvent")
    args = parser.parse_args()
    main(args.config_file, args.solvent_smiles)
    #python make_solvent.py ./config/make_sol.yaml "C1CCCCC1"
    #python make_solvent.py ./config/make_sol.yaml 'CC(=O)C'
    #python make_solvent.py ./config/make_sol.yaml 'CCO'
