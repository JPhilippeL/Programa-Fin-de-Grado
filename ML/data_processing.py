import os
from rdkit import Chem
import torch
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader

def mol_to_graph_data_obj(mol):
    bond_type_to_int = {
        Chem.rdchem.BondType.SINGLE: 0,
        Chem.rdchem.BondType.DOUBLE: 1,
        Chem.rdchem.BondType.TRIPLE: 2,
        Chem.rdchem.BondType.AROMATIC: 3
    }
    #Convierte una molécula RDKit a un objeto Data de PyTorch Geometric.
    
    atom_features = []
    for atom in mol.GetAtoms():
        atom_features.append([atom.GetAtomicNum()])
    x = torch.tensor(atom_features, dtype=torch.float)

    edge_index = []
    edge_attr = []
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()

        edge_index.append([i, j])
        edge_index.append([j, i])

        bond_type = bond.GetBondType()
        bond_type_idx = bond_type_to_int.get(bond_type, 0)  # default a 0

        edge_attr.append([bond_type_idx])
        edge_attr.append([bond_type_idx])

    edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
    edge_attr = torch.tensor(edge_attr, dtype=torch.float)

    data = Data(x=x, edge_index=edge_index, edge_attr=edge_attr)
    return data


def read_targets(targets_file):
    
    #Lee el archivo TXT de targets y devuelve un diccionario {mol_id: target}.
    
    target_dict = {}
    with open(targets_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            mol_id = parts[0]
            affinity = float(parts[1])
            target_dict[mol_id] = affinity
    return target_dict


def load_data_from_sdf(sdf_dir, target_dict):
    #Lee los archivos SDF de un directorio y los convierte a objetos Data con sus targets.

    data_list = []
    for filename in sorted(os.listdir(sdf_dir)):
        if filename.endswith('.sdf'):
            mol_id = filename.split('_')[0]
            if mol_id not in target_dict:
                print(f"Warning: No se encontró target para la molécula '{mol_id}', se omite.")
                continue

            affinity = target_dict[mol_id]
            filepath = os.path.join(sdf_dir, filename)
            suppl = Chem.SDMolSupplier(filepath, removeHs=False)
            mol = suppl[0]
            if mol is None:
                print(f"Warning: No se pudo leer la molécula del archivo '{filename}', se omite.")
                continue

            data = mol_to_graph_data_obj(mol)
            data.y = torch.tensor([affinity], dtype=torch.float)
            data_list.append(data)
    return data_list


def create_dataloader(data_list, batch_size=32, shuffle=True):
    return DataLoader(data_list, batch_size=batch_size, shuffle=shuffle)
