# main.py
import base64
import io
import os
from sdf_parser import parse_sdf
import graph_drawer
import os
import dash
from dash import html, dcc, Output, Input, State
import dash_cytoscape as cyto

import sdf_parser

"""
# === Configura el archivo SDF ===
#DEFAULT_FILENAME = "5RHD_ligand.sdf"
#FILE_PATH = os.path.join("data", DEFAULT_FILENAME)

# === Verifica si el archivo existe ===
#if not os.path.isfile(FILE_PATH):
#    raise FileNotFoundError(f"Archivo no encontrado: {FILE_PATH}")

# === Carga el grafo desde el archivo SDF ===
G = parse_sdf(FILE_PATH)
"""

# === Crea la app de Dash ===
app = dash.Dash(__name__)
app.title = "Visualizador de Moléculas SDF"

app.layout = html.Div([
    dcc.Upload(
        id='upload-data',
        children=html.Div(['Arrastra o haz click para seleccionar un archivo SDF']),
        style={
            'width': '100%', 'height': '60px', 'lineHeight': '60px',
            'borderWidth': '1px', 'borderStyle': 'dashed',
            'textAlign': 'center', 'margin': '10px',
            'cursor': 'pointer'
        },
        multiple=False,
        accept='.sdf'
    ),
    html.Div(id='message', children='Por favor, carga un archivo SDF para comenzar.'),
    cyto.Cytoscape(
        id='mol-graph',
        style={'width': '100%', 'height': '600px'},
        elements=[],
        stylesheet=graph_drawer.get_stylesheet(),
        layout={
            'name': 'cose',
            'randomize': False,
            'randomSeed': 1,
            'idealEdgeLength': 10,
            'nodeRepulsion': 400000,
            'edgeElasticity': 5,
            'gravity': 100,
            'numIter': 1000,
            'fit': True,
            'padding': 30
        }
    )
])

def parse_sdf_contents(contents):
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    with io.BytesIO(decoded) as f:
        G = parse_sdf(f)
        elements = graph_drawer.graph_to_cytoscape_elements(G)
    return elements

@app.callback(
    Output('mol-graph', 'elements'),
    Output('message', 'children'),
    Input('upload-data', 'contents')
)
def update_graph(contents):
    if contents is None:
        return [], 'Por favor, carga un archivo SDF para comenzar.'
    try:
        G = sdf_parser.parse_sdf_from_content(contents)
        elements = graph_drawer.graph_to_cytoscape_elements(G)
        return elements, ''
    except Exception as e:
        return [], f'Error al procesar el archivo: {str(e)}'

# === Iniciar el servidor ===
if __name__ == "__main__":
    app.run_server(debug=True)
