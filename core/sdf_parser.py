# sdf_parser.py
# Leer/Guardar archivos SDF y convertirlos a un grafo de NetworkX

import base64
from rdkit import Chem
from rdkit.Chem import AllChem
import networkx as nx

def parse_sdf(file_path):
    suppl = Chem.SDMolSupplier(file_path)
    mol = next((m for m in suppl if m is not None), None)

    if mol is None:
        raise ValueError("No se pudo leer una molécula válida desde el archivo SDF.")

    AllChem.Compute2DCoords(mol)  # asegura posiciones 2D
    graph = nx.Graph()

    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        pos = mol.GetConformer().GetAtomPosition(idx)
        graph.add_node(str(idx), 
               element=atom.GetSymbol(), 
               pos=(pos.x * 50, -pos.y * 50))  # escalar y reflejar eje Y si lo deseas

    for bond in mol.GetBonds():
        start = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        bond_type = bond.GetBondType()
        graph.add_edge(str(start), str(end), bond_type=str(bond_type))

    return graph

def parse_sdf_from_content(contents):
    # contents viene con un prefijo: "data:chemical/x-mdl-sdfile;base64,AAA..."
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)

    # rdkit espera un string con el contenido del archivo (mol block)
    mol_block = decoded.decode('utf-8')

    mol = Chem.MolFromMolBlock(mol_block, sanitize=True, removeHs=True)
    if mol is None:
        raise ValueError("No se pudo leer la molécula del contenido SDF.")
    
    graph = nx.Graph()

    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        pos = mol.GetConformer().GetAtomPosition(idx)
        graph.add_node(str(idx), 
                       element=atom.GetSymbol(), 
                       x=pos.x, y=pos.y, z=pos.z)

    for bond in mol.GetBonds():
        start = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        bond_type = bond.GetBondType()
        graph.add_edge(str(start), str(end), bond_type=str(bond_type))

    return graph
