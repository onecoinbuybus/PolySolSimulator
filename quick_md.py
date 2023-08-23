import os
import argparse
from radonpy.core import utils, poly
from radonpy.ff.gaff2_mod import GAFF2_mod
from radonpy.sim import qm
from polysol_md import MixMD, SamplingMD
from charge_utils import infer_charge_type


ff = GAFF2_mod()
ter_smiles = "*C"
omp_psi4 = 30
mpi = 30
omp = 30
gpu = 0
mem = 100000


def main(
    main_work_path,
    chain_length,
    in_monomer,
    in_esol,
    protocol="md1",
    min_dist_threshold=10,
):
    if not os.path.exists(main_work_path):
        os.makedirs(main_work_path)
        full_path = os.path.abspath(main_work_path)
        print(f"The working directory is: {full_path}")
    else:
        full_path = os.path.abspath(main_work_path)
        print(f"The working directory is already exists: {full_path}")

    work_dir = main_work_path
    polymer_mol = utils.JSONToMol(in_monomer)
    monomer_mol = utils.JSONToMol(in_monomer)
    charge_type = infer_charge_type(monomer_mol)
    if charge_type:
        print(f"The charge type set is: {charge_type}")
    else:
        print("Charge type could not be inferred.")

    ter = utils.mol_from_smiles(ter_smiles)
    qm.assign_charges(
        ter,
        charge=charge_type,
        opt=True,
        work_dir=work_dir,
        omp=omp_psi4,
        memory=mem,
        log_name="ter1",
    )
    polymer_mol = poly.polymerize_rw(monomer_mol, chain_length, tacticity="atactic")
    polymer_mol = poly.terminate_rw(polymer_mol, ter)
    result = ff.ff_assign(polymer_mol)

    equilibrated_solvent_mol = utils.JSONToMol(in_esol)
    mix_md = MixMD(equilibrated_solvent_mol, polymer_mol, work_dir, mpi, omp, gpu)
    equilibrated_mix = mix_md.exec(protocol, min_dist_threshold)

    sampling_md_instance = SamplingMD(work_dir, mpi, omp, gpu)
    sampling_result = sampling_md_instance.exec(equilibrated_mix)
    print("finish")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulation Script")
    parser.add_argument("main_work_path", type=str, help="Path to main work directory")
    parser.add_argument("chain_length", type=int, help="Chain length for the polymer")
    parser.add_argument(
        "-in_monomer", type=str, required=True, help="Path to monomer JSON file"
    )
    parser.add_argument(
        "-in_esol",
        type=str,
        required=True,
        help="Path to equilibrated solvent JSON file",
    )
    parser.add_argument(
        "-protocol", type=str, default="md1", help="Simulation protocol"
    )
    parser.add_argument(
        "-min_dist_threshold", type=int, default=10, help="Minimum distance threshold"
    )

    args = parser.parse_args()
    main(
        args.main_work_path,
        args.chain_length,
        args.in_monomer,
        args.in_esol,
        args.protocol,
        args.min_dist_threshold
    )
    # python quick_md.py poly30_ethanol 30 -in_monomer ps.json -in_esol ethanol.json



