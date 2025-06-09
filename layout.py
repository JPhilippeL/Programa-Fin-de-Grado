# layout.py
from dash import html, dcc
import dash_cytoscape as cyto
import graph_drawer

homelayout = html.Div([
    dcc.Upload(
        id='upload-data',
        children=html.Div(['Nuevo'], style={'display': 'inline-block'}),
        style={
            'position': 'absolute',  # para posicionarlo relativo a la ventana
            'top': '10px',           # distancia desde arriba
            'left': '10px',          # distancia desde la izquierda
            'width': '100px',        # ancho más pequeño, parecido a un botón
            'height': '40px',
            'lineHeight': '40px',
            'borderWidth': '1px',
            'borderStyle': 'solid',
            'borderRadius': '5px',   # bordes redondeados como botón
            'textAlign': 'center',
            'backgroundColor': '#007bff',  # color azul típico botón
            'color': 'white',
            'cursor': 'pointer',
            'userSelect': 'none',
            'fontWeight': 'bold',
            'zIndex': 1000  # para asegurarse de que esté encima de otros elementos
        },
        multiple=False,
        accept='.sdf'
    ),
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
            'zIndex': 1
        },
        userPanningEnabled=False,  # deshabilita mover la cámara
        userZoomingEnabled=True,   # permite hacer zoom
        minZoom=0.5,               # zoom mínimo (ajusta según quieras)
        maxZoom=2                  # zoom máximo (ajusta según quieras)
    )
])
