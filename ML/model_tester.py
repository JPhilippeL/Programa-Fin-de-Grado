#model_tester.py

import torch
from rdkit import Chem
from torch_geometric.data import Data
from ML.model_trainer import create_model
from ML.data_processing import mol_to_graph_data_obj
import os
import glob
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error
from math import sqrt
from ML.data_processing import read_targets, load_data_from_sdf
import logging
logger = logging.getLogger(__name__)
logging.getLogger('matplotlib.font_manager').setLevel(logging.WARNING)
logging.getLogger('matplotlib').setLevel(logging.WARNING)

def cargar_y_predecir(checkpoint_path, sdf_path):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # Cargar checkpoint
    checkpoint = torch.load(checkpoint_path, map_location=device)

     # Recuperar metadatos
    model_type = checkpoint['model_type']
    input_dim = checkpoint['input_dim']
    edge_dim = checkpoint['edge_dim']
    hidden_dim = checkpoint.get('hidden_dim', 64)
    num_layers = checkpoint.get('num_layers', 3)
    target_name = checkpoint.get('target_name', 'target')

    # Crear modelo con los parámetros guardados
    model = create_model(
        model_type,
        input_dim=input_dim,
        edge_dim=edge_dim,
        hidden_dim=hidden_dim,
        num_layers=num_layers
    )
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

    # Predecir
    with torch.no_grad():
        batch = torch.zeros(data.num_nodes, dtype=torch.long, device=data.x.device)  # todos nodos del mismo grafo
        out = model(data.x, data.edge_index, data.edge_attr, batch)
        pred = out.squeeze().item()

    return pred, target_name

def test_model_on_directory(checkpoint_path, sdf_dir, targets_file):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # Crear carpeta de resultados si no existe
    resultados_dir = "Resultados"
    os.makedirs(resultados_dir, exist_ok=True)

    # Cargar checkpoint
    checkpoint = torch.load(checkpoint_path, map_location=device)
    model = create_model(
        checkpoint['model_type'],
        input_dim=checkpoint['input_dim'],
        edge_dim=checkpoint['edge_dim'],
        hidden_dim=checkpoint.get('hidden_dim', 64),
        num_layers=checkpoint.get('num_layers', 3)
    )
    model.load_state_dict(checkpoint['model_state_dict'])
    model.to(device)
    model.eval()

    # Leer datos
    target_dict = read_targets(targets_file)
    data_list = load_data_from_sdf(sdf_dir, target_dict)

    y_true, y_pred, filenames = [], [], []

    for data in data_list:
        data = data.to(device)
        batch = torch.zeros(data.num_nodes, dtype=torch.long, device=device)

        with torch.no_grad():
            out = model(data.x, data.edge_index, data.edge_attr, batch)
            pred = out.squeeze().item()

        y_pred.append(pred)
        y_true.append(data.y.item())
        filenames.append(data.name if hasattr(data, 'name') else 'unknown')

     # Nombre base de archivos
    model_filename = os.path.basename(checkpoint_path)
    model_name_no_ext = os.path.splitext(model_filename)[0]
    folder_name = os.path.basename(sdf_dir.rstrip(os.sep))

    # Archivo de predicciones
    output_predictions_path = os.path.join(
        resultados_dir,
        f"predicciones_{model_name_no_ext}_{folder_name}.txt"
    )

    with open(output_predictions_path, 'w') as f:
        for fname, pred in zip(filenames, y_pred):
            f.write(f"{fname} {pred:.4f}\n")

    # RMSE
    rmse = sqrt(mean_squared_error(y_true, y_pred))
    logger.info(f"RMSE: {rmse:.4f}")

    # Scatter plot
    plt.figure(figsize=(6, 6))
    plt.scatter(y_true, y_pred, alpha=0.7)
    plt.plot([min(y_true), max(y_true)], [min(y_true), max(y_true)], color='red', linestyle='--')
    plt.xlabel("Valor real")
    plt.ylabel("Predicción")
    plt.title(f"Scatter Plot - {model_name_no_ext} - {folder_name}")
    plt.grid(True)
    plt.tight_layout()

    # Guardar imagen
    plot_filename = os.path.join(
        resultados_dir,
        f"scatter_plot_{model_name_no_ext}_{folder_name}.png"
    )
    plt.savefig(plot_filename)
    plt.close()

    logger.info(f"Scatter plot guardado en: {plot_filename}")
    logger.info(f"Predicciones guardadas en: {output_predictions_path}")
