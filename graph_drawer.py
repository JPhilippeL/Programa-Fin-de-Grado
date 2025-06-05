import matplotlib.pyplot as plt
from matplotlib import collections as mc
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout

def draw_graph_with_pygraphviz(G, prog='neato'):
    """
    Dibuja el grafo G usando PyGraphviz para calcular el layout.
    `prog` puede ser: 'dot', 'neato', 'fdp', 'circo', 'twopi'...
    """
    pos = graphviz_layout(G, prog=prog)

    draw_bonds(G, pos)  # Asumiendo que tu función personalizada usa pos igual que nx

    nx.draw_networkx_nodes(G, pos, node_color='lightgreen', node_size=225)
    nx.draw_networkx_labels(G, pos, labels={n: G.nodes[n].get('element', str(n)) for n in G.nodes})

    plt.axis("equal")
    plt.title(f"Grafo molecular con layout '{prog}' de Graphviz")
    plt.show()

def draw_bonds(G, pos):
    lines = []
    styles = []

    for u, v in G.edges():
        bond_type = G[u][v]['bond_type']
        x1, y1 = pos[u]
        x2, y2 = pos[v]

        if bond_type == "SINGLE":
            lines.append([(x1, y1), (x2, y2)])
            styles.append(('solid', 1.5))

        elif bond_type == "DOUBLE":
            offset = 4
            dx, dy = y2 - y1, x1 - x2
            norm = (dx**2 + dy**2)**0.5 or 1e-6
            dx, dy = dx / norm * offset, dy / norm * offset
            lines.append([(x1 + dx, y1 + dy), (x2 + dx, y2 + dy)])
            lines.append([(x1 - dx, y1 - dy), (x2 - dx, y2 - dy)])
            styles.extend([('solid', 1.5)] * 2)

        elif bond_type == "TRIPLE":
            offset = 4
            dx, dy = y2 - y1, x1 - x2
            norm = (dx**2 + dy**2)**0.5 or 1e-6
            dx, dy = dx / norm * offset, dy / norm * offset
            lines.append([(x1, y1), (x2, y2)])
            lines.append([(x1 + dx, y1 + dy), (x2 + dx, y2 + dy)])
            lines.append([(x1 - dx, y1 - dy), (x2 - dx, y2 - dy)])
            styles.extend([('solid', 1.5)] * 3)

        else:  # Por ejemplo, aromático o desconocido
            lines.append([(x1, y1), (x2, y2)])
            styles.append(('dotted', 1.0))

    for line, (style, width) in zip(lines, styles):
        lc = mc.LineCollection([line], linestyles=style, linewidths=width, colors='black')
        plt.gca().add_collection(lc)
