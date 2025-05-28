# sdf_parser.py

from rdkit import Chem
import networkx as nx

def parse_sdf(file_path):
    suppl = Chem.SDMolSupplier(file_path)
    mol = next((m for m in suppl if m is not None), None)

    if mol is None:
        raise ValueError("No se pudo leer una molécula válida desde el archivo SDF.")

    graph = nx.Graph()

    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        pos = mol.GetConformer().GetAtomPosition(idx)
        graph.add_node(idx, 
                       element=atom.GetSymbol(), 
                       x=pos.x, y=pos.y, z=pos.z)

    for bond in mol.GetBonds():
        start = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        bond_type = bond.GetBondType()
        graph.add_edge(start, end, bond_type=str(bond_type))

    return graph
