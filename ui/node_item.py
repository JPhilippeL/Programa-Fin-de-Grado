from PySide6.QtWidgets import QGraphicsEllipseItem, QGraphicsTextItem, QGraphicsItem, QMenu
from PySide6.QtGui import QBrush, QColor, QFont
from PySide6.QtCore import Signal, QObject
from ui.utils import ATOM_COLORS, ATOM_TEXT_COLORS, ATOM_COLORS_DEFAULT, ATOM_TEXT_COLORS_DEFAULT

NODE_RADIUS = 20

class NodeItem(QGraphicsEllipseItem, QObject):
    def __init__(self, x, y, radius, element, node_id):
        super().__init__(-radius, -radius, 2 * radius, 2 * radius)
        
         # Estilo según el elemento
        color = QColor(ATOM_COLORS.get(element, ATOM_COLORS_DEFAULT))
        text_color = QColor(ATOM_TEXT_COLORS.get(element, ATOM_TEXT_COLORS_DEFAULT))

        self.node_id = node_id
        self.setBrush(QBrush(QColor(color)))
        self.setFlag(QGraphicsItem.ItemIsMovable)
        self.setFlag(QGraphicsItem.ItemSendsGeometryChanges)
        self.setZValue(1)

        # Texto del elemento (como "C", "O", etc.)
        self.label = QGraphicsTextItem(element, self)
        font = QFont("Arial", 12)
        font.setBold(True)
        self.label.setFont(font)

        self.label.setDefaultTextColor(text_color)
        self.label.setZValue(2)
        # Calcular centro del nodo y ajustar el texto
        text_rect = self.label.boundingRect()
        self.label.setPos(-text_rect.width() / 2, -text_rect.height() / 2)

        # Lista de conexiones (enlaces)
        self.edges = []

        # Posición inicial
        self.setPos(x, y)

    def itemChange(self, change, value):
        if change == QGraphicsItem.ItemPositionChange:
            for edge in self.edges:
                edge.update_position()
        return super().itemChange(change, value)
    
    modify_requested = Signal(object)  # self
    delete_requested = Signal(object)  # self
    add_edge_requested = Signal(object)  # self

    def contextMenuEvent(self, event):
        menu = QMenu()
        modify_action = menu.addAction("Modificar nodo")
        delete_action = menu.addAction("Eliminar nodo")
        add_edge_action = menu.addAction("Añadir enlace")

        action = menu.exec(event.screenPos())

        if action == modify_action:
            self.modify_requested.emit(self)
        elif action == delete_action:
            self.delete_requested.emit(self)
        elif action == add_edge_action:
            self.add_edge_requested.emit(self)