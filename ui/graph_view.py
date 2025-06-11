from PySide6.QtWidgets import QGraphicsView, QGraphicsScene, QGraphicsEllipseItem, QGraphicsLineItem, QGraphicsTextItem, QGraphicsItem, QGraphicsItemGroup
from PySide6.QtGui import QPen, QBrush, QColor, QPainter, QFont
from PySide6.QtCore import Qt
from ui.utils import ATOM_COLORS, ATOM_TEXT_COLORS, ATOM_COLORS_DEFAULT, ATOM_TEXT_COLORS_DEFAULT, BACKGROUND_COLOR

NODE_RADIUS = 20

class NodeItem(QGraphicsEllipseItem):
    def __init__(self, x, y, radius, element):
        super().__init__(-radius, -radius, 2 * radius, 2 * radius)
        
         # Estilo según el elemento
        color = QColor(ATOM_COLORS.get(element, ATOM_COLORS_DEFAULT))
        text_color = QColor(ATOM_TEXT_COLORS.get(element, ATOM_TEXT_COLORS_DEFAULT))

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


class EdgeItem(QGraphicsItemGroup):
    def __init__(self, source_node, target_node, bond_type):
        super().__init__()
        self.source = source_node
        self.target = target_node
        self.bond_type = bond_type.upper()

        self.lines = []
        self._create_lines()
        self.update_position()

        self.source.edges.append(self)
        self.target.edges.append(self)

    def _create_lines(self):
        # Determina desplazamientos según el tipo de enlace
        if self.bond_type == 'DOUBLE':
            offsets = [-2, 2]
        elif self.bond_type == 'TRIPLE':
            offsets = [-4, 0, 4]
        else:
            offsets = [0]  # SINGLE o AROMATIC

        for offset in offsets:
            line = QGraphicsLineItem()
            pen = QPen(Qt.white, 2)

            if self.bond_type == 'AROMATIC':
                pen.setStyle(Qt.DashLine)

            line.setPen(pen)
            line.setZValue(0)
            line.offset = offset  # desplazamiento perpendicular
            self.addToGroup(line)
            self.lines.append(line)

    def update_position(self):
        src = self.source.sceneBoundingRect().center()
        tgt = self.target.sceneBoundingRect().center()

        dx = tgt.x() - src.x()
        dy = tgt.y() - src.y()
        length = (dx ** 2 + dy ** 2) ** 0.5
        if length == 0:
            return

        # Vector perpendicular normalizado
        norm_x = -dy / length
        norm_y = dx / length

        for line in self.lines:
            offset = line.offset
            ox = norm_x * offset
            oy = norm_y * offset

            line.setLine(
                src.x() + ox, src.y() + oy,
                tgt.x() + ox, tgt.y() + oy
            )



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


class MoleculeGraphView(QGraphicsView):
    def __init__(self, graph):
        super().__init__()
        self.setScene(MoleculeGraphScene(graph))
        self.setRenderHint(QPainter.Antialiasing)
        self.setDragMode(QGraphicsView.ScrollHandDrag)
        self.setViewportUpdateMode(QGraphicsView.FullViewportUpdate)

