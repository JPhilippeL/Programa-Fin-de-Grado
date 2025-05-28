# main.py

import argparse
import os
from sdf_parser import parse_sdf
import matplotlib.pyplot as plt
import networkx as nx

def draw_graph(G):
    pos = {node: (G.nodes[node]['x'], G.nodes[node]['y']) for node in G.nodes}
    labels = {node: G.nodes[node]['element'] for node in G.nodes}
    edge_labels = {(u, v): G[u][v]['bond_type'] for u, v in G.edges}

    nx.draw(G, pos, labels=labels, with_labels=True, node_size=600, node_color="lightgreen")
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)
    plt.title("Molécula desde SDF con RDKit")
    plt.axis("equal")
    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Visualizar molécula desde archivo SDF.")
    parser.add_argument("filename", help="Nombre del archivo SDF (sin la ruta)")
    args = parser.parse_args()

    file_path = os.path.join("data", args.filename)

    if not os.path.isfile(file_path):
        print(f"Archivo no encontrado: {file_path}")
        return

    G = parse_sdf(file_path)
    draw_graph(G)

if __name__ == "__main__":
    main()
