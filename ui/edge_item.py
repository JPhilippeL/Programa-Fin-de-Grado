from PySide6.QtWidgets import QGraphicsLineItem, QGraphicsItemGroup
from PySide6.QtGui import QPen
from PySide6.QtCore import Qt

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