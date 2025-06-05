# main.py

import argparse
import os
from sdf_parser import parse_sdf
import graph_drawer
import matplotlib.pyplot as plt
import networkx as nx


# === Convertir a formato Cytoscape ===
def nx_to_cytoscape(G):
    nodes = [
        {
            'data': {'id': str(n), 'label': G.nodes[n].get('element', str(n))},
        }
        for n in G.nodes
    ]
    edges = [
        {
            'data': {'source': str(u), 'target': str(v)},
        }
        for u, v in G.edges
    ]
    return nodes + edges

def main():
    parser = argparse.ArgumentParser(description="Visualizar molécula desde archivo SDF.")
    parser.add_argument("filename", nargs='?', default="5RGV_ligand.sdf", help="Nombre del archivo SDF")
    args = parser.parse_args()

    file_path = os.path.join("data", args.filename)

    if not os.path.isfile(file_path):
        print(f"Archivo no encontrado: {file_path}")
        return

    G = parse_sdf(file_path)
    graph_drawer.draw_graph_with_pygraphviz(G, prog='neato')
    #draw_graph(G)

if __name__ == "__main__":
    main()
