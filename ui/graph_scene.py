from PySide6.QtWidgets import QGraphicsScene, QInputDialog, QMessageBox
from PySide6.QtGui import QColor
from ui.utils import BACKGROUND_COLOR
from ui.edge_item import EdgeItem
from ui.node_item import NodeItem
from core.graph_manager import GraphManager

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
            node_id = data.get("idx", node_id)
            node = NodeItem(x, y, 20, element, node_id)

            node.modify_requested.connect(self.on_modify_node)
            node.delete_requested.connect(self.on_delete_node)
            node.add_edge_requested.connect(self.on_add_edge)

            self.addItem(node)
            self.node_items[node_id] = node

        # Crear enlaces
        for source_id, target_id, data in self.graph.edges(data=True):
            bond_type = data.get('bond_type').upper()
            source_node = self.node_items[source_id]
            target_node = self.node_items[target_id]
            edge = EdgeItem(source_node, target_node, bond_type)
            self.addItem(edge)

    def on_modify_node(self, node_item):
        node_id = self._get_node_id_from_item(node_item)
        new_element, ok = QInputDialog.getText(None, "Modificar nodo", "Nuevo símbolo químico:")
        if ok and new_element:
            GraphManager.modify_node(self.graph, node_id, new_element)
            node_item.label.setPlainText(new_element.upper())

    def on_delete_node(self, node_item):
        node_id = self._get_node_id_from_item(node_item)
        confirm = QMessageBox.question(None, "Confirmar", "¿Eliminar este nodo?")
        if confirm == QMessageBox.Yes:
            GraphManager.delete_node(self.graph, node_id)
            for edge in node_item.edges[:]:
                self.removeItem(edge)
            self.removeItem(node_item)
            del self.node_items[node_id]

    def on_add_edge(self, node_item):
        # Aquí puedes implementar selección de segundo nodo
        print(f"Añadir enlace desde {self._get_node_id_from_item(node_item)}")

    def _get_node_id_from_item(self, item):
        # Asume que puedes mapear el item a su id
        for node_id, node_item in self.node_items.items():
            if node_item == item:
                return node_id