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
输入polymer smiles以及电荷种类，以下面步骤创建高分子单体:  
1. **构象搜索**: 使用RDKit搜索分子的合适构象。
2. **分子预处理**: 对分子进行预处理，建立trimer
3. **电荷赋值**: 使用psi4对分子电荷（gasteiger, RESP, ESP, Mulliken, Lowdin）进行计算，提取trimer中间部分电荷，建立monomer

可以通过以下命令创建monomer:

```bash
python mol_gen.py "*C(C*)c1ccccc1" -dir "test_run" -out "ps.json"
```
下面是可用的命令以及描述:
下面是可用的命令行参数及其描述：
- `smiles`: polymer smiles字符串
- `-dir`: 工作目录名
- `-out`: 输出文件名，用于保存生成的分子JSON文件。
- `-charge`: 电荷类型，默认为"RESP"。

### 2. Make solvent
创建溶剂盒子并进行预平衡。可以指定溶剂分子，电荷以及密度。此过程包括以下步骤:  
1. **构象搜索**: 使用RDKit搜索溶剂分子的合适构象。
2. **电荷赋值**: 使用psi4对溶剂分子的电荷进行计算。
3. **创建溶剂单元**: 根据指定的分子数量和密度，创建一个无定型或简单的溶剂单元。
溶剂盒子可用以下命令创建：
```bash
python make_solvent.py ./config/make_sol.yaml "C1CCCCC1"
```
下面是命令的具体描述
- `config_file`: YAML配置文件路径，其中包括OMP线程数、MPI进程数、GPU使用量、内存、溶剂分子数量、密度和盒子类型等参数。
- `solvent_smiles`: 溶剂分子的smiles
### 3. Execute Quick Molecular Dynamics 

这个脚本用于制作高分子-溶剂混合盒子并进行模拟以及抽样计算。可以通过以下命令执行快速MD模拟:

```bash
python quick_md.py poly30_ethanol 30 -in_monomer ps.json -in_esol ethanol.json
```

下面是可用的命令行参数及其描述：
- `main_work_path`: 主工作目录的路径。
- `chain_length`: 聚合物的链长。
- `-in_monomer`: 单体JSON文件的路径。
- `-in_esol`: 平衡溶剂JSON文件的路径。
- `-protocol`: 模拟协议，md1: 8步的压缩与解压缩。md2: 保持压力与温度进行模拟。默认为"md1"。  
- `-min_dist_threshold`: 初始高分子与溶剂间的最小距离阈值，默认为10。

## Details of Simulation

### Step 1. Monomer preparation  
使用三个重复链块的中间部分获得了部分电荷(RESP)并且创建了指定高分子单体的Mol文件。RESP使用psi4进行计算。单体使用RDKit的etkdgv2进行构象搜索和优化。之后使用Random Walk创建高分子单链用于模拟。

### Step 2. Equilibrate simulation of solvent  
随机生成初始密度为 0.05g/cm³的溶剂盒子，里面放入1000个溶剂分子。使用NPT模拟令纯溶剂盒以2fs的时间步长均匀平衡2 ns。温度固定为300K，压力为1 atm。

### Step 3. Creation of initial mixture cell  
将生成的单链polymer放入溶剂盒的中心，删除了与高分子重叠以及距离过近的溶剂分子。

### Step 4. Simulation of polymer-solvent system  
固定高分子的端点原子，使用NPT模拟了2,000,000步的混合物系统，步长为5fs。温度固定为300K，这里使用了并且比较了两个方案:
  - 方案1: 将压强升高到100atm之后解压到1atm。
  - 方案2: 保持压强为1直到模拟结束。

### Step 5. Sampling  
松开端点原子，以300K，压力为1 atm的条件以同样的时间步长模拟500,000步的NPT并进行各种性质的计算。

## Result

### Reference
Kangjie LYU, Yanqiu PENG, Li XIAO, Juntao LU, Lin ZHUANG. Atomistic Understanding of the Peculiar Dissolution Behavior of Alkaline Polymer Electrolytes in Alcohol/Water Mixed Solvents. _Acta Phys. -Chim. Sin._, 2019, Vol. 35, Issue (4): 378-384. doi: [10.3866/PKU.WHXB201805031](https://doi.org/10.3866/PKU.WHXB201805031)

Yoshihiro Hayashi, Junichiro Shiomi, Junko Morikawa & Ryo Yoshida. RadonPy: automated physical property calculation using all-atom classical molecular dynamics simulations for polymer informatics. _npj Computational Materials_, volume 8, Article number: 222 (2022). [Cite this article]([https://doi.org/10.1038/s41524-022-00906-4)https://doi.org/10.1038/s41524-022-00906-4]）


