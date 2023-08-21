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
## Advanced Usage
```python

```
## Details of simulation
Step 1. Monomer preparation
使用三个重复链块的中间部分获得了部分电荷(RESP)并且创建了指定高分子单体的Mol文件。RESP使用psi4进行计算。单体使用RDKit的etkdgv2进行构象搜索和优化。之后使用Random Walk创建高分子单链用于模拟。

Step 2. Equilibrate simulation of solvent
随机生成初始密度为 0.05g/cm³的溶剂盒子，里面放入1000个溶剂分子。使用NPT模拟令纯溶剂盒以2fs的时间步长均匀平衡2 ns。温度固定为300K，压力为1 atm。

Step 3. Creation of initial mixture cell
将生成的单链polymer放入溶剂盒的中心，删除了与高分子重叠以及距离过近的溶剂分子。

Step 4. Simulation of polymer-solvent system
固定高分子的端点原子，使用NPT模拟了2。000，000步的混合物系统，步长为5fs。温度固定为300K，这里使用了并且比较了两个方案:
方案1: 将压强升高到100 atm之后解压到1atm。
方案2: 保持压强为1直到模拟结束。

Step 5. Sampling
松开端点原子，以300K，压力为1 atm的条件以同样的时间步长模拟500,000步的NPT并进行各种性质的计算。


