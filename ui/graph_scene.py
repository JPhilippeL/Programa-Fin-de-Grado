from PySide6.QtWidgets import QGraphicsScene, QInputDialog, QMessageBox, QGraphicsLineItem
from PySide6.QtGui import QColor, QPen, Qt, QTransform
from PySide6.QtCore import Signal, QObject
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

        self.edge_mode = False
        self.edge_source_node = None
        self.temp_edge_line = None

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
            node.add_edge_requested.connect(self.start_temporary_edge)

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
            node_item.update_element(new_element)

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
            
    def start_temporary_edge(self, source_node):
        self.edge_mode = True
        self.edge_source_node = source_node
        self.temp_edge_line = QGraphicsLineItem()
        self.temp_edge_line.setPen(QPen(Qt.white, 2, Qt.SolidLine))
        self.addItem(self.temp_edge_line)

    def cancel_temporary_edge(self):
        if self.temp_edge_line:
            self.removeItem(self.temp_edge_line)
            self.temp_edge_line = None
        self.edge_mode = False
        self.edge_source_node = None

    def mouseMoveEvent(self, event):
        if self.edge_mode and self.temp_edge_line:
            start = self.edge_source_node.sceneBoundingRect().center()
            end = event.scenePos()
            self.temp_edge_line.setLine(start.x(), start.y(), end.x(), end.y())
        super().mouseMoveEvent(event)

    def mousePressEvent(self, event):
        if self.edge_mode:
            if event.button() == Qt.RightButton:
                self.cancel_temporary_edge()
                return

            if event.button() == Qt.LeftButton:
                clicked_item = self.itemAt(event.scenePos(), QTransform())
                # Subir la jerarquía hasta encontrar un NodeItem
                while clicked_item and not isinstance(clicked_item, NodeItem):
                    clicked_item = clicked_item.parentItem()

                if isinstance(clicked_item, NodeItem) and clicked_item != self.edge_source_node:
                    GraphManager.add_edge(
                        self.graph,
                        self._get_node_id_from_item(self.edge_source_node),
                        self._get_node_id_from_item(clicked_item),
                        bond_type="SINGLE"
                    )
                    newedge = EdgeItem(self.edge_source_node, clicked_item, "SINGLE")
                    self.addItem(newedge)
                    self.cancel_temporary_edge()
                    return


        super().mousePressEvent(event)

    def keyPressEvent(self, event):
        if self.edge_mode and event.key() == Qt.Key_Escape:
            self.cancel_temporary_edge()
            return
        super().keyPressEvent(event)