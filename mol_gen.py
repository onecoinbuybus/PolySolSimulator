import argparse
from radonpy.core import utils
from radonpy.ff.gaff2_mod import GAFF2_mod
from radonpy.sim import qm
from polysol_md import make_work_dir
from charge_utils import AtomIndexing, set_charges, get_selected_charges_list

def main(smiles, work_dir_name, output_file_name, charge_type="RESP"):
    ff = GAFF2_mod()
    omp_psi4 = 10
    mpi = 10
    omp = 10
    mem = 100000

    main_work_path = make_work_dir(work_dir_name)
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
        log_name="monomer1",
        etkdg_ver=2,
    )

    mol_update = utils.deepcopy(mol)

    atom_indexing = AtomIndexing(mol)
    sorted_result_dict, processed_mol, homopoly_rw = atom_indexing.process_mol()

    qm.assign_charges(
        homopoly_rw,
        charge=charge_type,
        opt=False,
        work_dir=work_dir,
        omp=omp_psi4,
        memory=mem,
        log_name="monomer1",
    )

    selected_indexes = list(sorted_result_dict.values())
    selected_charges = get_selected_charges_list(
        homopoly_rw, selected_indexes, charge=charge_type
    )

    mol_update = set_charges(mol_update, selected_charges, charge_type=charge_type)
    utils.MolToJSON(mol_update, output_file_name)
    print("finish mol gen, file name:", output_file_name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run main program")
    parser.add_argument("smiles", help="SMILES string")
    parser.add_argument("-dir", dest="work_dir_name", help="Work directory name", required=True)
    parser.add_argument("-out", dest="output_file_name", help="Output file name", required=True)
    parser.add_argument("-charge", dest="charge", help="Charge type", default="RESP")  # 定义-charge作为可选参数
    args = parser.parse_args()
    main(args.smiles, args.work_dir_name, args.output_file_name, args.charge)

    # python mol_gen.py "*C(C*)c1ccccc1" -dir "test_run" -out "ps.json" -charge "gasteiger"
    # python mol_gen.py "*C(C*)c1ccccc1" -dir "test_run" -out "ps.json" (RESP)



