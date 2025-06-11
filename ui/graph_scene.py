from PySide6.QtWidgets import QGraphicsScene
from PySide6.QtGui import QColor
from ui.utils import BACKGROUND_COLOR
from ui.edge_item import EdgeItem
from ui.node_item import NodeItem

class MoleculeGraphScene(QGraphicsScene):
    def __init__(self, graph):
        super().__init__()
        self.setBackgroundBrush(QColor(BACKGROUND_COLOR))
        self.graph = graph
        self._draw_graph()

    def _draw_graph(self):
        self.node_items = {}
        
        # Crear nodos
        for node_id, data in self.graph.nodes(data=True):
            x, y = data["pos"]
            element = data["element"]
            node = NodeItem(x, y, 20, element)
            self.addItem(node)
            self.node_items[node_id] = node

        # Crear enlaces
        for source_id, target_id, data in self.graph.edges(data=True):
            bond_type = data.get('bond_type').upper()
            source_node = self.node_items[source_id]
            target_node = self.node_items[target_id]
            edge = EdgeItem(source_node, target_node, bond_type)
            self.addItem(edge)