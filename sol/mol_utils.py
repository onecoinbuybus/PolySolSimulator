from collections import Counter
from collections import defaultdict
from rdkit import Chem


def assign_original_indices_to_atoms(mol):
    """
    Set 'Mid_unit' and 'Ori_idx' properties for atoms in a given molecule.

    Args:
        mol: RDKit Mol object
    """
    ori_idx_info = []
    try:
        for i in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(i)
            # 使用 SetProp 方法为原子设置一个名字
            atom.SetProp("Mid_unit", "mid_atom")
            original_index = str(atom.GetIdx())
            atom.SetProp("Ori_idx", original_index)
            ori_idx_info.append(original_index)

        return ori_idx_info
    except:
        print("fail to set idx on atom")


# 得到trimer里面的原子序号信息
def extract_atom_information(molecule):
    """
    Function to get atom information from the given molecule.

    Args:
        molecule: RDKit Mol object

    Returns:
        atom_info: A list of dictionaries containing atom information
    """
    atom_info = []
    for i in range(molecule.GetNumAtoms()):
        atom = molecule.GetAtomWithIdx(i)
        try:
            atom.GetProp("Mid_unit")
            original_index = atom.GetProp("Ori_idx")
            symbol = atom.GetSymbol()
            atom_index = atom.GetIdx()
            atom_info.append(
                {"atom": atom_index, "original_index": original_index, "symbol": symbol}
            )
        except:
            pass
    return atom_info


def group_atoms_by_original_index(atom_info):
    atoms_by_original_index = defaultdict(list)
    for atom in atom_info:
        atoms_by_original_index[atom["original_index"]].append(atom["atom"])
    return atoms_by_original_index


def find_occurrences_by_original_index(atoms_by_original_index):
    single_occurrence_original_indices = [
        original_index
        for original_index, atoms in atoms_by_original_index.items()
        if len(atoms) == 1
    ]
    second_occurrence_by_original_index = {
        original_index: atoms[1]
        for original_index, atoms in atoms_by_original_index.items()
        if len(atoms) > 1
    }
    return single_occurrence_original_indices, second_occurrence_by_original_index


def construct_sorted_result_by_original_index(
    single_occurrence_original_indices,
    second_occurrence_by_original_index,
    atoms_by_original_index,
):
    result_dict = {}
    for original_index in single_occurrence_original_indices:
        result_dict[original_index] = atoms_by_original_index[original_index][0]
    for original_index, atom_index in second_occurrence_by_original_index.items():
        result_dict[original_index] = atom_index
    sorted_result_by_original_index = {
        str(k): v for k, v in sorted(result_dict.items(), key=lambda item: int(item[0]))
    }
    return sorted_result_by_original_index


def get_charges(mol, charge="gasteiger"):
    charges_by_index = {}

    if charge == "gasteiger":
        Chem.rdPartialCharges.ComputeGasteigerCharges(mol)
        for i, atom in enumerate(mol.GetAtoms()):
            charges_by_index[i] = {
                "index": i,
                "AtomicCharge": float(atom.GetProp("_GasteigerCharge")),
            }

    elif charge == "RESP":
        for i, atom in enumerate(mol.GetAtoms()):
            charges_by_index[i] = {
                "index": i,
                "RESP": atom.GetDoubleProp("RESP"),
                "ESP": atom.GetDoubleProp("ESP"),
                "AtomicCharge": atom.GetDoubleProp("RESP"),
            }

    elif charge == "ESP":
        for i, atom in enumerate(mol.GetAtoms()):
            charges_by_index[i] = {
                "index": i,
                "RESP": atom.GetDoubleProp("RESP"),
                "ESP": atom.GetDoubleProp("ESP"),
                "AtomicCharge": atom.GetDoubleProp("ESP"),
            }

    elif charge == "Mulliken":
        for i, atom in enumerate(mol.GetAtoms()):
            charges_by_index[i] = {
                "index": i,
                "MullikenCharge": atom.GetDoubleProp("MullikenCharge"),
                "AtomicCharge": atom.GetDoubleProp("MullikenCharge"),
            }

    elif charge == "Lowdin":
        for i, atom in enumerate(mol.GetAtoms()):
            charges_by_index[i] = {
                "index": i,
                "LowdinCharge": atom.GetDoubleProp("LowdinCharge"),
                "AtomicCharge": atom.GetDoubleProp("LowdinCharge"),
            }

    else:
        print(f"{charge} is not implemented.")
        return None

    return charges_by_index


def set_charges(mol, charges_list, charge_type="gasteiger"):
    """
    Assigns atomic charges to RDKit Mol object based on the provided charges list.

    Args:
        mol: RDKit Mol object
        charges_list: A list containing dictionaries with atomic charges information
        charge_type: The type of charge to be assigned (default: 'gasteiger')

    Returns:
        RDKit Mol object with assigned charges
    """

    if charge_type == "gasteiger":
        for i in range(len(charges_list)):
            mol.GetAtomWithIdx(i).SetDoubleProp(
                "AtomicCharge", float(charges_list[i]["AtomicCharge"])
            )

    elif charge_type == "RESP":
        for i in range(len(charges_list)):
            mol.GetAtomWithIdx(i).SetDoubleProp("RESP", float(charges_list[i]["RESP"]))
            mol.GetAtomWithIdx(i).SetDoubleProp("ESP", float(charges_list[i]["ESP"]))
            mol.GetAtomWithIdx(i).SetDoubleProp(
                "AtomicCharge", float(charges_list[i]["AtomicCharge"])
            )

    elif charge_type == "ESP":
        for i in range(len(charges_list)):
            mol.GetAtomWithIdx(i).SetDoubleProp("RESP", charges_list[i]["RESP"])
            mol.GetAtomWithIdx(i).SetDoubleProp("ESP", charges_list[i]["ESP"])
            mol.GetAtomWithIdx(i).SetDoubleProp(
                "AtomicCharge", charges_list[i]["AtomicCharge"]
            )

    elif charge_type == "Mulliken":
        for i in range(len(charges_list)):
            mol.GetAtomWithIdx(i).SetDoubleProp(
                "MullikenCharge", charges_list[i]["MullikenCharge"]
            )
            mol.GetAtomWithIdx(i).SetDoubleProp(
                "AtomicCharge", charges_list[i]["AtomicCharge"]
            )

    elif charge_type == "Lowdin":
        for i in range(len(charges_list)):
            mol.GetAtomWithIdx(i).SetDoubleProp(
                "LowdinCharge", charges_list[i]["LowdinCharge"]
            )
            mol.GetAtomWithIdx(i).SetDoubleProp(
                "AtomicCharge", charges_list[i]["AtomicCharge"]
            )

    else:
        print(f"{charge_type} is not implemented.")
        return None

    return mol


def get_selected_charges_list(mol, selected_indexes, charge="gasteiger"):
    selected_charges_list = []

    # 获取所有原子的电荷信息
    charges_dict = get_charges(mol, charge)

    # 仅选择 selected_indexes 中指定的原子的电荷信息
    for atom_index in selected_indexes:
        selected_charges_list.append(charges_dict[atom_index])

    return selected_charges_list
