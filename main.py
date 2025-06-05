# main.py
import argparse
import os
from sdf_parser import parse_sdf
import graph_drawer
import matplotlib.pyplot as plt
import networkx as nx
import os
import dash
from dash import html, Output, Input
import dash_cytoscape as cyto

# === Configura el archivo SDF ===
DEFAULT_FILENAME = "5RHD_ligand.sdf"
FILE_PATH = os.path.join("data", DEFAULT_FILENAME)

# === Verifica si el archivo existe ===
if not os.path.isfile(FILE_PATH):
    raise FileNotFoundError(f"Archivo no encontrado: {FILE_PATH}")

# === Carga el grafo desde el archivo SDF ===
G = parse_sdf(FILE_PATH)

# === Crea la app de Dash ===
app = dash.Dash(__name__)
app.title = "Visualizador de Moléculas SDF"

app.layout = html.Div([
    cyto.Cytoscape(
        id='mol-graph',
        style={'width': '100%', 'height': '600px'},
        elements= graph_drawer.graph_to_cytoscape_elements(G),
        stylesheet= graph_drawer.get_stylesheet(),
        layout = {
            'name': 'cose',
            'randomize': False,         # Usa el orden actual de los nodos como punto de partida
            'randomSeed': 1,           # Fija la semilla
            'idealEdgeLength': 10,     # Distancia "ideal" entre nodos conectados
            'nodeRepulsion': 400000,    # Fuerza de repulsión entre nodos
            'edgeElasticity': 5,      # Qué tan elástica es la arista (fuerza)
            'gravity': 100,              # Atracción general hacia el centro
            'numIter': 1000,            # Cuántas iteraciones ejecutar (más = más estable)
            'fit': True,                # Ajusta el grafo al contenedor
            'padding': 30               # Margen desde los bordes
        }

    ),
    html.Div(id='output')
])

# === Iniciar el servidor ===
if __name__ == "__main__":
    app.run_server(debug=True)
