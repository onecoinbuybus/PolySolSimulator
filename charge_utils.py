from radonpy.core import poly
from collections import defaultdict
from rdkit import Chem

class AtomIndexing:
    def __init__(self, mol):
        self.mol = mol

    def set_atom_idx(self):
        ori_idx_info = []
        try:
            for i in range(self.mol.GetNumAtoms()):
                atom = self.mol.GetAtomWithIdx(i)
                atom.SetProp("Mid_unit", "mid_atom")
                original_index = str(atom.GetIdx())
                atom.SetProp("Ori_idx", original_index)
                ori_idx_info.append(original_index)

            return ori_idx_info
        except:
            print("fail to set idx on atom")

    def get_atom_info(self):
        atoms_info = []
        for idx in range(self.mol.GetNumAtoms()):
            atom = self.mol.GetAtomWithIdx(idx)
            try:
                original_index = atom.GetProp("Ori_idx")
                symbol = atom.GetSymbol()
                atom_index = atom.GetIdx()
                atoms_info.append(
                    {
                        "index": atom_index,
                        "original_index": original_index,
                        "symbol": symbol,
                    }
                )
            except:
                pass
        return atoms_info

    def get_sorted_result_dict(self):
        atom_info = self.get_atom_info()
        original_index_atoms = defaultdict(list)

        for atom in atom_info:
            original_index_atoms[atom["original_index"]].append(atom["index"])

        single_occurrence_original_indexes = [
            original_index
            for original_index, atoms in original_index_atoms.items()
            if len(atoms) == 1
        ]
        second_occurrence_atoms = {
            original_index: atoms[1]
            for original_index, atoms in original_index_atoms.items()
            if len(atoms) > 1
        }

        result_dict = {}
        for original_index in single_occurrence_original_indexes:
            result_dict[original_index] = original_index_atoms[original_index][0]

        for original_index, atom_index in second_occurrence_atoms.items():
            result_dict[original_index] = atom_index

        sorted_result_dict = {
            str(k): v
            for k, v in sorted(result_dict.items(), key=lambda item: int(item[0]))
        }

        return sorted_result_dict

    def process_mol(self):
        _ = self.set_atom_idx()
        homopoly = poly.polymerize_rw(self.mol, 3)
        indexing = AtomIndexing(homopoly)
        sorted_dict = indexing.get_sorted_result_dict()
        return sorted_dict, self.mol, homopoly


def get_charges(mol, charge="gasteiger"):
    charges_dict = {}

    if charge == "gasteiger":
        Chem.rdPartialCharges.ComputeGasteigerCharges(mol)
        for i, atom in enumerate(mol.GetAtoms()):
            charges_dict[i] = {
                "index": i,
                "AtomicCharge": float(atom.GetProp("_GasteigerCharge")),
            }

    elif charge == "RESP":
        for i, atom in enumerate(mol.GetAtoms()):
            charges_dict[i] = {
                "index": i,
                "RESP": atom.GetDoubleProp("RESP"),
                "ESP": atom.GetDoubleProp("ESP"),
                "AtomicCharge": atom.GetDoubleProp("RESP"),
            }

    elif charge == "ESP":
        for i, atom in enumerate(mol.GetAtoms()):
            charges_dict[i] = {
                "index": i,
                "RESP": atom.GetDoubleProp("RESP"),
                "ESP": atom.GetDoubleProp("ESP"),
                "AtomicCharge": atom.GetDoubleProp("ESP"),
            }

    elif charge == "Mulliken":
        for i, atom in enumerate(mol.GetAtoms()):
            charges_dict[i] = {
                "index": i,
                "MullikenCharge": atom.GetDoubleProp("MullikenCharge"),
                "AtomicCharge": atom.GetDoubleProp("MullikenCharge"),
            }

    elif charge == "Lowdin":
        for i, atom in enumerate(mol.GetAtoms()):
            charges_dict[i] = {
                "index": i,
                "LowdinCharge": atom.GetDoubleProp("LowdinCharge"),
                "AtomicCharge": atom.GetDoubleProp("LowdinCharge"),
            }

    else:
        print(f"{charge} is not implemented.")
        return None

    return charges_dict


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

    charges_dict = get_charges(mol, charge)

    for atom_index in selected_indexes:
        selected_charges_list.append(charges_dict[atom_index])

    return selected_charges_list

