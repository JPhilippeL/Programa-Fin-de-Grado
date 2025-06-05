#graph_drawer.py
# Convertir un grafo de NetworkX a formato Cytoscape para visualización
import utils

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

def graph_to_cytoscape_elements(G):
    elements = []
    for node, data in G.nodes(data=True):
        element = data['element']
        color = utils.ATOM_COLORS.get(element, utils.ATOM_COLORS_DEFAULT)
        text_color = utils.ATOM_TEXT_COLORS.get(element, utils.ATOM_TEXT_COLORS_DEFAULT)  # negro por defecto si no existe
        elements.append({
            'data': {
                'id': str(node),
                'label': element,
                'color': color,
                'text_color': text_color
            }
        })
    for u, v, data in G.edges(data=True):
        elements.append({
            'data': {
                'source': str(u),
                'target': str(v),
                'bond_type': data.get('bond_type', 'SINGLE')
            }
        })
    return elements


def get_stylesheet():
    stylesheet = [
    {
        'selector': 'node',
        'style': {
            'content': 'data(label)',
            'background-color':'data(color)',   # verde más oscuro y vibrante
            'color': 'data(text_color)',                 # color del texto claro
            'font-size': '14px',
            'text-valign': 'center',
            'text-halign': 'center',
            'width': 30,
            'height': 30
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
    return stylesheet
