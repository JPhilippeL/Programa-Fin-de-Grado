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

stylesheet = [
    {
        'selector': 'node',
        'style': {
            'content': 'data(label)',
            'background-color': '#4CAF50',   # verde más oscuro y vibrante
            'color': '#eee',                 # color del texto claro
            'font-size': '14px',
            'text-valign': 'center',
            'text-halign': 'center',
            'width': 30,
            'height': 30,
            'border-width': 2,
            'border-color': '#2e7d32'       # borde verde oscuro
        }
    },
    {
        'selector': 'edge',
        'style': {
            'curve-style': 'bezier',
            'line-color': '#E9ECF5',            # gris medio para enlaces por defecto
            'width': 2,
            'target-arrow-shape': 'none',
        }
    },
    # SINGLE
    {
        'selector': '[bond_type = "SINGLE"]',
        'style': {
            'line-style': 'solid',
            'width': 2
        }
    },
    # DOUBLE
    {
        'selector': '[bond_type = "DOUBLE"]',
        'style': {
            'line-style': 'solid',
            'width': 5
        }
    },
    # TRIPLE
    {
        'selector': '[bond_type = "TRIPLE"]',
        'style': {
            'line-style': 'solid',
            'width': 7
        }
    },
    # AROMATIC or unknown
    {
        'selector': '[bond_type = "AROMATIC"]',
        'style': {
            'line-style': 'dashed',
            'width': 2
        }
    },
    # Texto de nodos en tema oscuro
    {
        'selector': 'node > text',
        'style': {
            'color': '#eee'
        }
    }
]

# === Configura el archivo SDF ===
DEFAULT_FILENAME = "5RGV_ligand.sdf"
FILE_PATH = os.path.join("data", DEFAULT_FILENAME)

# === Verifica si el archivo existe ===
if not os.path.isfile(FILE_PATH):
    raise FileNotFoundError(f"Archivo no encontrado: {FILE_PATH}")

# === Carga el grafo desde el archivo SDF ===
G = parse_sdf(FILE_PATH)

# === Convierte a formato Cytoscape ===
def nx_to_cytoscape(G):
    nodes = [
        {'data': {'id': str(n), 'label': G.nodes[n].get('element', str(n))}}
        for n in G.nodes
    ]
    edges = [
        {
            'data': {
                'source': str(u),
                'target': str(v),
                'bond_type': str(G.edges[u, v].get('bond_type', '')).upper(),  # FORZAR MAYÚSCULAS
                'label': G.edges[u, v].get('bond_type', '')
            }
        }
        for u, v in G.edges
    ]
    return nodes + edges

# === Crea la app de Dash ===
app = dash.Dash(__name__)
app.title = "Visualizador de Moléculas SDF"

app.layout = html.Div([
    cyto.Cytoscape(
        id='mol-graph',
        layout={'name': 'cose', 'randomSeed': 3},  # Puedes probar 'cose', 'circle', 'breadthfirst', etc.
        style={'width': '100%', 'height': '600px'},
        elements=nx_to_cytoscape(G),
        stylesheet=stylesheet,
    ),
    html.Div(id='output')
])

# === Iniciar el servidor ===
if __name__ == "__main__":
    app.run_server(debug=True)
