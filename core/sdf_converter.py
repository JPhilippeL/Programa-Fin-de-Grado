# sdf_parser.py
# Leer/Guardar archivos SDF y convertirlos a un grafo de NetworkX

from rdkit import Chem
from rdkit.Chem import AllChem
import networkx as nx

def parse_sdf(file_path):
    suppl = Chem.SDMolSupplier(file_path, removeHs=False)
    mol = next((m for m in suppl if m is not None), None)

    if mol is None:
        raise ValueError("No se pudo leer una molécula válida desde el archivo SDF.")
    
    AllChem.Compute2DCoords(mol)
    graph = nx.Graph()

    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        pos = mol.GetConformer().GetAtomPosition(idx)
        graph.add_node(str(idx), 
               element=atom.GetSymbol(), 
               pos=(pos.x * 50, -pos.y * 50))   # Multiplicamos por 50 para escalar las coordenadas
    for bond in mol.GetBonds():
        start = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        bond_type = bond.GetBondType()
        graph.add_edge(str(start), str(end), bond_type=str(bond_type))

    return graph

def graph_to_mol(graph):
    mol = Chem.RWMol()

    # Mapear ID de nodos a nuevos índices de átomos
    node_to_idx = {}

    for node_id in sorted(graph.nodes, key=int):
        element = graph.nodes[node_id]["element"]
        atom = Chem.Atom(element)
        idx = mol.AddAtom(atom)
        node_to_idx[node_id] = idx

    for u, v, data in graph.edges(data=True):
        bond_type_str = data.get("bond_type", "SINGLE").upper()
        bond_type = getattr(Chem.rdchem.BondType, bond_type_str, Chem.rdchem.BondType.SINGLE)
        mol.AddBond(node_to_idx[u], node_to_idx[v], bond_type)

    mol = mol.GetMol()
    AllChem.Compute2DCoords(mol)

    return mol

def save_graph_as_sdf(graph, file_path):
    mol = graph_to_mol(graph)
    try:
        Chem.SanitizeMol(mol)
        writer = Chem.SDWriter(file_path)
        writer.write(mol)
        writer.close()
    except Exception as e:
        raise RuntimeError(f"Error al guardar la molécula: {str(e)}")