import MDAnalysis as mda
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from rdkit import Chem
import os

def plot_rg(rg_file_path, save_path='rg.png'):
    df_sam = pd.read_csv(rg_file_path, skiprows=4, delim_whitespace=True, header=None, names=['id', 'c_gyr1'])
    rgs = df_sam[df_sam['id']==1]['c_gyr1']    
    median_value = np.median(rgs)
    plt.plot(rgs, label='rg')
    plt.axhline(y=median_value, color='r', linestyle='--', label='Median')
    plt.xlabel('Step')
    plt.ylabel('rg')
    plt.legend()
    plt.savefig(save_path)
    absolute_path = os.path.abspath(save_path)
    print(f"Rg plot has been saved at: {absolute_path}")
    plt.show()


def get_carbon_atoms_by_residue_names(poly, residue_names=['TU0', 'TU1']):
    # 存储与指定残基名称匹配的碳原子的索引
    carbon_atom_indices = []

    # 遍历mol中的所有原子
    for idx, atom in enumerate(poly.GetAtoms()):
        if atom.GetSymbol() == 'C':  # 检查是否为碳原子
            residue_info = atom.GetPDBResidueInfo()
            if residue_info is not None:
                res_name = residue_info.GetResidueName()
                if res_name in residue_names:
                    carbon_atom_indices.append(idx)

    return carbon_atom_indices


def plot_distance(dump_file_path, h_idx, t_idx, save_path='end_to_end_distance.png'):   
    universe = mda.Universe(dump_file_path, format="LAMMPSDUMP")
    distances = []
    
    for ts in universe.trajectory:
        distance_vector = universe.atoms[h_idx].position - universe.atoms[t_idx].position
        distance = np.linalg.norm(distance_vector)
        distances.append(distance)
        
    median_value = np.median(distances)
    
    plt.plot(distances, label='end-to-end distance')
    plt.axhline(y=median_value, color='r', linestyle='--', label='Median')
    plt.xlabel('Step (*1000)')
    plt.ylabel('distance')
    plt.legend()
    plt.savefig(save_path)  # 保存图像为PNG文件
    plt.savefig(save_path)
    absolute_path = os.path.abspath(save_path)
    print(f"End-to-end distance plot has been saved at: {absolute_path}")
    plt.show()
    plt.show()



def display_png(image_path):
    img = mpimg.imread(image_path)
    plt.imshow(img)
    plt.axis('off')  
    plt.show()

