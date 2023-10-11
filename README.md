# PolySolSimulator
Automatic Tool for All-Atom Molecular Dynamics Simulations of Polymer-Solvent Systems Based on Radonpy

## Introduction
PolySolSimulator is an automated library for running equilibration molecular dynamics (EMD) simulations of polymer-solvent mixed systems. By inputting the SMILES of the polymer and solvent, along with parameters such as polymer length, the library automatically creates a mixed cell and conducts simulations. It can calculate the radius of gyration (Rg) and end-to-end distance of the polymer within the solvent.

## Requirements
- Python 3.8
- LAMMPS
- rdkit
- psi4
- mdtraj
- scipy

## Simple Usage
```bash
git clone https://github.com/onecoinbuybus/PolySolSimulator/
cd PolySolSimulator
python run_md.py
```
Firstly a file named "test_run" will be created.The solvent cell of actone will be automatically generated, and then combined with polystyrene to form a mixed cell for EMD simulation  

## Advanced Usage
### 1. Create monomer
Input the polymer smiles and charge type, then follow the steps below to create a polymer monomer:
1.**Conformational search:** Use RDKit to search for appropriate conformations of the molecule.
2.**Molecular preprocessing:** Preprocess the molecule and construct a trimer.
3.**Charge assignment:** Use psi4 to calculate molecular charges (gasteiger, RESP, ESP, Mulliken, Lowdin) and extract the central part of the charge in the trimer to construct the monomer.


Create monomer using following command:

```bash
python mol_gen.py "*C(C*)c1ccccc1" -dir "test_run" -out "ps.json"
```


### 2. Make solvent
Create a solvent box and pre-equilibrate it. You can specify the solvent molecule, charge, and density. The process includes the following steps:

1.Conformational search: Use RDKit to search for appropriate conformations of the solvent molecule.
2.Charge assignment: Use psi4 to calculate charges for the solvent molecule.
3.Construct solvent unit: Based on the specified molecular quantity and density, create an amorphous or simple solvent unit.
4.The solvent box can be created with the following command:
    
```bash
python make_solvent.py ./config/make_sol.yaml "C1CCCCC1"
```

### 3. Execute Quick Molecular Dynamics 

This script is used to create a polymer-solvent mixed box and carry out simulations as well as sampling calculations. You can perform a quick MD simulation using the following command:

```bash
python quick_md.py poly30_ethanol 30 -in_monomer ps.json -in_esol ethanol.json
```

## Available command-line parameters and their descriptions:

- **main_work_path**: Path to the main working directory.
- **chain_length**: Chain length of the polymer.
- **-in_monomer**: Path to the monomer JSON file.
- **-in_esol**: Path to the equilibrated solvent JSON file.
- **-protocol**: Simulation protocol. 
  - `md1`: 8-step compression and decompression. 
  - `md2`: Simulation maintaining pressure and temperature. 
  - Default is "md1".
- **-min_dist_threshold**: Initial minimum distance threshold between the polymer and solvent, default is 10.


## Details of Simulation

### Step 1. Monomer preparation  
The partial charges (RESP) were obtained using the middle section of three repeating chain blocks, and a Mol file of the specified polymer monomer was created. RESP calculations were performed with psi4. The monomer underwent conformational search and optimization using RDKit's etkdgv2. A polymer single chain for simulation was then created using Random Walk.

### Step 2. Equilibrate simulation of solvent  
A solvent box with an initial density of 0.05g/cm³ was randomly generated, containing 1,000 solvent molecules. An NPT simulation was conducted to equilibrate the pure solvent box for 2 ns with a timestep of 2fs. The temperature was maintained at 300K, and the pressure at 1 atm.

### Step 3. Creation of initial mixture cell  
The generated single-chain polymer was placed at the center of the solvent box, and solvent molecules overlapping with the polymer or too close to it were removed.

### Step 4. Simulation of polymer-solvent system  
With the end atoms of the polymer fixed, an NPT simulation of the mixture system was performed for 2,000,000 steps with a timestep of 5fs. The temperature was maintained at 300K. Two scenarios were used and compared:
  - Scenario 1: The pressure was increased to 100 atm and then decompressed to 1 atm.
  - Scenario 2: The pressure was maintained at 1 atm until the end of the simulation.

### Step 5. Sampling  
Releasing the endpoint atoms, an NPT simulation was conducted for 500,000 steps with the same timestep under conditions of 300K and a pressure of 1 atm, and various properties were calculated.


## Result

### Reference
Kangjie LYU, Yanqiu PENG, Li XIAO, Juntao LU, Lin ZHUANG. Atomistic Understanding of the Peculiar Dissolution Behavior of Alkaline Polymer Electrolytes in Alcohol/Water Mixed Solvents. _Acta Phys. -Chim. Sin._, 2019, Vol. 35, Issue (4): 378-384. doi: [10.3866/PKU.WHXB201805031](https://doi.org/10.3866/PKU.WHXB201805031)

Yoshihiro Hayashi, Junichiro Shiomi, Junko Morikawa & Ryo Yoshida. RadonPy: automated physical property calculation using all-atom classical molecular dynamics simulations for polymer informatics. _npj Computational Materials_, volume 8, Article number: 222 (2022). [Cite this article]([https://doi.org/10.1038/s41524-022-00906-4)https://doi.org/10.1038/s41524-022-00906-4]）


