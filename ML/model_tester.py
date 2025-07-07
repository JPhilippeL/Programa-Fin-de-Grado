#model_tester.py

import torch
from rdkit import Chem
from torch_geometric.data import Data
from ML.model_trainer import create_model
from ML.data_processing import mol_to_graph_data_obj

def cargar_y_predecir(checkpoint_path, sdf_path):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # Cargar checkpoint
    checkpoint = torch.load(checkpoint_path, map_location=device)

    # Recuperar metadatos
    model_type = checkpoint['model_type']
    input_dim = checkpoint['input_dim']
    edge_dim = checkpoint['edge_dim']
    target_name = checkpoint.get('target_name', 'target')

    # Crear modelo con los parámetros guardados
    model = create_model(model_type, input_dim, edge_dim)
    model.load_state_dict(checkpoint['model_state_dict'])
    model.to(device)
    model.eval()

    # Leer molécula del SDF
    suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
    mol = suppl[0]
    if mol is None:
        raise ValueError(f"No se pudo leer la molécula de {sdf_path}")

    # Convertir a PyG Data
    data = mol_to_graph_data_obj(mol)
    data = data.to(device)

    # Predecir (suponiendo que el modelo devuelve un tensor con la predicción)
    with torch.no_grad():
        batch = torch.zeros(data.num_nodes, dtype=torch.long, device=data.x.device)  # todos nodos del mismo grafo
        out = model(data.x, data.edge_index, data.edge_attr, batch)
        pred = out.squeeze().item()

    return pred, target_name
