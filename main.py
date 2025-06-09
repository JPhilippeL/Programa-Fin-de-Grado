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
from layout import homelayout  # Importamos el layout separado


import sdf_parser

app = dash.Dash(__name__)
app.title = "Visualizador de Moléculas SDF"

app.layout = homelayout

def parse_sdf_contents(contents):
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    with io.BytesIO(decoded) as f:
        G = parse_sdf(f)
        elements = graph_drawer.graph_to_cytoscape_elements(G)
    return elements

@app.callback(
    Output('mol-graph', 'elements'),
    Input('upload-data', 'contents')
)
def update_graph(contents):
    if contents is None:
        return []
    try:
        G = sdf_parser.parse_sdf_from_content(contents)
        elements = graph_drawer.graph_to_cytoscape_elements(G)
        return elements
    except Exception as e:
        return [], f'Error al procesar el archivo: {str(e)}'

# === Iniciar el servidor ===
if __name__ == "__main__":
    app.run_server(debug=True)
