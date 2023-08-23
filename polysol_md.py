import os
import cell_utils
from radonpy.core import utils, poly
from radonpy.sim.lammps import LAMMPS
from radonpy.sim.lammps import MolFromLAMMPSdata
from radonpy.sim.md import MD
from datetime import datetime
import numpy as np


def make_work_dir(main_file_name):
    work_dir_path = os.path.abspath(main_file_name)
    try:
        os.makedirs(os.path.join(work_dir_path, "sol_md"))
        os.makedirs(os.path.join(work_dir_path, "mix_md"))
        os.makedirs(os.path.join(work_dir_path, "sampling_md"))
        print(f"Directories created successfully under {work_dir_path}/")
        return work_dir_path
    except:
        print("Directory already exists.")
        return work_dir_path


class Poly_Sol_EMD:
    def __init__(
        self, solvent_cell, polymer_mol, work_dir, ff, omp_psi4, mpi, omp, mem, gpu
    ):
        self.solvent_cell = solvent_cell
        self.polymer_mol = polymer_mol
        self.work_dir = work_dir
        self.ff = ff
        self.omp_psi4 = omp_psi4
        self.mpi = mpi
        self.omp = omp
        self.mem = mem
        self.gpu = gpu

    def sol_emd(
        self,
        sol_work_dir=None,
        cutoff_in=12.0,
        cutoff_out=16.0,
        p_start=1,
        p_stop=1,
        t_start=300,
        t_stop=300,
        step=1000000,
        time_step=0.002,
        shake=True,
    ):
        if sol_work_dir is None:
            sol_work_dir = os.path.join(self.work_dir, "sol_md")
            os.makedirs(sol_work_dir, exist_ok=True)

        md = MD()
        md.cutoff_in = cutoff_in
        md.cutoff_out = cutoff_out
        md.dat_file = "sol_cell_ini.data"
        md.dump_file = "sol_npt.dump"
        md.xtc_file = "sol_npt.xtc"
        md.rst = True
        md.write_data = "sol_npt_last.data"
        md.add_min(min_style="cg")
        md.add_md(
            "npt",
            step,
            time_step=time_step,
            t_start=t_start,
            t_stop=t_stop,
            p_start=p_start,
            p_stop=p_stop,
            shake=shake,
        )
        return md

    def poly_sol_emd1(
        self,
        cutoff_in=12.0,
        cutoff_out=16.0,
        press_ratio=None,
        start_press=1,
        max_press=1000,
        step_weight=None,
        step=2000000,
        time_step=0.005,
        t_start=300,
        t_stop=300,
        p_dump=500,
    ):
        md = MD()
        md.cutoff_in = cutoff_in
        md.cutoff_out = cutoff_out
        md.dat_file = "mix_cell_ini.data"
        md.dump_file = "poly_sol_npt.dump"
        md.xtc_file = "poly_sol.xtc"
        md.rst = True
        md.write_data = "poly_sol_last.data"
        md.add.append("group fixed_atoms id 1")  # 包括你想要固定的原子
        md.add.append("group rest_atoms subtract all fixed_atoms")  # 其他所有原子
        md.add.append("fix fix1 fixed_atoms setforce 0.0 0.0 0.0")
        if press_ratio is None:
            press_ratio = [0.02, 0.60, 1.00, 0.50, 0.10, 0.01]

        # NPT compression/decompression equilibration protocol, 8 step
        press_list = (
            [start_press] + list(np.array(press_ratio) * max_press) + [start_press]
        )
        # step list
        if step_weight is None:
            step_weight = [0.05, 0.1, 0.3, 0.2, 0.1, 0.05, 0.2]
        step_list = [int(i * step) for i in step_weight]
        for s, p in zip(step_list, press_list):
            md.add_md(
                "npt",
                group="rest_atoms",
                step=s,
                time_step=time_step,
                t_start=t_start,
                t_stop=t_stop,
                p_start=p,
                p_stop=p,
                p_dump=p_dump,
                add=["neigh_modify delay 0 every 1 check no"],
                shake=False,
            )
        return md

    def poly_sol_emd2(
        self,
        cutoff_in=12.0,
        cutoff_out=16.0,
        step=2000000,
        time_step=0.005,
        p_start=1,
        p_stop=1,
        p_dump=500,
        t_start=300,
        t_stop=300,
    ):
        md = MD()
        md.cutoff_in = cutoff_in
        md.cutoff_out = cutoff_out
        md.dat_file = "mix_cell_ini.data"
        md.dump_file = "poly_sol_npt.dump"
        md.xtc_file = "poly_sol.xtc"
        md.rst = True
        md.write_data = "poly_sol_last.data"
        md.add.append("group fixed_atoms id 1")  # 包括你想要固定的原子
        md.add.append("group rest_atoms subtract all fixed_atoms")  # 其他所有原子
        md.add.append("fix fix1 fixed_atoms setforce 0.0 0.0 0.0")
        md.add_md(
            "npt",
            group="rest_atoms",
            step=step,
            time_step=time_step,
            t_start=t_start,
            t_stop=t_stop,
            p_start=p_start,
            p_stop=p_stop,
            p_dump=p_dump,
            shake=False,
        )
        return md

    def sampling(
        self,
        data_file=None,
        cutoff_in=12.0,
        cutoff_out=16.0,
        step=500000,
        time_step=0.005,
        p_start=1,
        p_stop=1,
        p_dump=500,
        t_start=300,
        t_stop=300,
    ):
        mix_md_dir = os.path.join(self.work_dir, "mix_md")
        if data_file is None:
            data_file = os.path.join(mix_md_dir, "poly_sol_last.data")
        else:
            data_file=data_file
        md = MD()
        md.cutoff_in = cutoff_in
        md.cutoff_out = cutoff_out
        md.dat_file = data_file
        md.dump_file = "sampling.dump"
        md.xtc_file = "sampling.xtc"
        md.rst = True
        md.write_data = "sampling_last.data"
        md.add_md(
            "npt",
            group="all",
            step=step,
            time_step=time_step,
            t_start=t_start,
            t_stop=t_stop,
            p_start=p_start,
            p_stop=p_stop,
            p_dump=p_dump,
            shake=True,
        )
        md.wf[-1].add_rg()
        return md

    def exec(self, protocol, min_dist_threshold=5):
        sol_work_dir = os.path.join(self.work_dir, "sol_md")
        mix_work_dir = os.path.join(self.work_dir, "mix_md")
        sampling_work_dir = os.path.join(self.work_dir, "sampling")

        # solvent equilibration
        print("start solvent equilibration")

        solvent_cell_c = utils.deepcopy(self.solvent_cell)
        lmp_sol = LAMMPS(work_dir=sol_work_dir)
        lmp_sol.make_dat(solvent_cell_c, confId=0, file_name="sol_cell_ini.data")
        sol_md = self.sol_emd()
        lmp_sol.make_input(sol_md)
        equilibrated_sol = lmp_sol.run(
            sol_md, mol=solvent_cell_c, omp=self.omp, mpi=self.mpi, gpu=self.gpu
        )
        utils.MolToJSON(
            equilibrated_sol, os.path.join(sol_work_dir, "equilibrated_solvent.json")
        )

        # mix md
        print("start mixture equilibration")
        ini_mix_cell = cell_utils.create_mixture_cell(
            polymer_mol=self.polymer_mol,
            solvent_cell=equilibrated_sol,
            min_dist_threshold=min_dist_threshold,
        )
        ini_mix_cell_c = utils.deepcopy(ini_mix_cell)
        lmp_mix = LAMMPS(work_dir=mix_work_dir)
        lmp_mix.make_dat(ini_mix_cell_c, confId=0, file_name="mix_cell_ini.data")
        if protocol == "md1":
            md_mix_emd = self.poly_sol_emd1()
        elif protocol == "md2":
            md_mix_emd = self.poly_sol_emd2()
        lmp_mix.make_input(md_mix_emd)
        equilibrated_mix = lmp_mix.run(
            md_mix_emd, mol=ini_mix_cell_c, omp=self.omp, mpi=self.mpi, gpu=self.gpu
        )
        utils.MolToJSON(
            equilibrated_mix, os.path.join(mix_work_dir, "equilibrated_mix.json")
        )

        # sampling
        print("start sampling")
        equilibrated_mix_c = utils.deepcopy(equilibrated_mix)
        lmp_sampling = LAMMPS(work_dir=sampling_work_dir)
        lmp_sampling.make_dat(
            equilibrated_mix_c, confId=0, file_name="equilibrated_cell_ini.data"
        )
        sampling_md = self.sampling()
        lmp_sampling.make_input(sampling_md)
        sampling_mix = lmp_sampling.run(
            sampling_md,
            mol=equilibrated_mix_c,
            omp=self.omp,
            mpi=self.mpi,
            gpu=self.gpu,
        )
        utils.MolToJSON(
            sampling_mix, os.path.join(sampling_work_dir, "sampling_mix.json")
        )

        return sampling_mix


class SolventEquilibration(Poly_Sol_EMD):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def exec(self):
        sol_work_dir = os.path.join(self.work_dir, "sol_md")
        os.makedirs(sol_work_dir, exist_ok=True)
        start_time = datetime.now()
        # solvent equilibration
        print("start solvent equilibration", "time:", start_time)
        solvent_cell_c = utils.deepcopy(self.solvent_cell)
        lmp_sol = LAMMPS(work_dir=sol_work_dir)
        lmp_sol.make_dat(solvent_cell_c, confId=0, file_name="sol_cell_ini.data")
        sol_md = self.sol_emd()
        lmp_sol.make_input(sol_md)
        equilibrated_sol = lmp_sol.run(
            sol_md, mol=solvent_cell_c, omp=self.omp, mpi=self.mpi, gpu=self.gpu
        )
        print(
            "finish npt solvent equilibration, time cost:", datetime.now() - start_time
        )
        utils.MolToJSON(
            equilibrated_sol, os.path.join(sol_work_dir, "equilibrated_solvent.json")
        )


class MixMD(Poly_Sol_EMD):
    def __init__(self, solvent_cell, polymer_mol, work_dir, mpi, omp, gpu):
        self.solvent_cell = solvent_cell
        self.polymer_mol = polymer_mol
        self.work_dir = work_dir
        self.mpi = mpi
        self.omp = omp
        self.gpu = gpu

    def exec(self, protocol, min_dist_threshold=5):
        mix_work_dir = os.path.join(self.work_dir, "mix_md")
        os.makedirs(mix_work_dir, exist_ok=True)
        start_time = datetime.now()
        print("start mixture equilibration", "time:", start_time)
        ini_mix_cell = cell_utils.create_mixture_cell(
            polymer_mol=self.polymer_mol,
            solvent_cell=self.solvent_cell,
            min_dist_threshold=min_dist_threshold,
        )
        ini_mix_cell_c = utils.deepcopy(ini_mix_cell)
        lmp_mix = LAMMPS(work_dir=mix_work_dir)
        lmp_mix.make_dat(ini_mix_cell_c, confId=0, file_name="mix_cell_ini.data")
        if protocol == "md1":
            md_mix_emd = self.poly_sol_emd1()
        elif protocol == "md2":
            md_mix_emd = self.poly_sol_emd2()
        lmp_mix.make_input(md_mix_emd)
        equilibrated_mix = lmp_mix.run(
            md_mix_emd, mol=ini_mix_cell_c, omp=self.omp, mpi=self.mpi, gpu=self.gpu
        )
        print(
            "finish npt mixture equilibration, time cost:", datetime.now() - start_time
        )
        utils.MolToJSON(
            equilibrated_mix, os.path.join(mix_work_dir, "equilibrated_mix.json")
        )
        return equilibrated_mix


class SamplingMD(Poly_Sol_EMD):
    def __init__(self, work_dir, mpi, omp, gpu):
        self.work_dir = os.path.abspath(work_dir)
        self.mpi = mpi
        self.omp = omp
        self.gpu = gpu
        
    def exec(self, equilibrated_mix_mol):
        sampling_work_dir = os.path.join(self.work_dir, "sampling_md")
        os.makedirs(sampling_work_dir, exist_ok=True)
        start_time = datetime.now()
        print("start sampling", "time:", start_time)
        equilibrated_mix_c = utils.deepcopy(equilibrated_mix_mol)
        lmp_sampling = LAMMPS(work_dir=sampling_work_dir)
        lmp_sampling.make_dat(
            equilibrated_mix_c, confId=0, file_name="equilibrated_mix_cell_ini.data"
        )
        data_file_path = os.path.join(sampling_work_dir, "equilibrated_mix_cell_ini.data")
        sampling_md = self.sampling(data_file=data_file_path)
        lmp_sampling.make_input(sampling_md)
        sampling_mix = lmp_sampling.run(
            sampling_md,
            mol=equilibrated_mix_c,
            omp=self.omp,
            mpi=self.mpi,
            gpu=self.gpu,
        )
        print("finish sampling, time cost:", datetime.now() - start_time)
        return sampling_mix
