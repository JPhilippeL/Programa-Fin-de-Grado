from PySide6.QtWidgets import QMenuBar, QFileDialog
from PySide6.QtGui import QAction
from core.sdf_parser import parse_sdf
from ui.graph_view import MoleculeGraphView

class MenuBar(QMenuBar):
    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent  # Referencia a MainWindow

        archivo_menu = self.addMenu("Archivo")

        cargar_action = QAction("Cargar", self)
        cargar_action.triggered.connect(self.cargar_archivo)
        archivo_menu.addAction(cargar_action)

        guardar_action = QAction("Guardar", self)
        archivo_menu.addAction(guardar_action)  # sin funcionalidad por ahora

    def cargar_archivo(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self.parent,
            "Seleccionar archivo SDF",
            "",
            "Archivos SDF (*.sdf);;Todos los archivos (*)"
        )
        if file_path:
            graph = parse_sdf(file_path)

            if self.parent.graph_view:
                self.parent.centralWidget().layout().removeWidget(self.parent.graph_view)
                self.parent.graph_view.deleteLater()

            self.parent.graph_view = MoleculeGraphView(graph)

            if not self.parent.centralWidget():
                from PySide6.QtWidgets import QWidget, QVBoxLayout
                central = QWidget()
                layout = QVBoxLayout(central)
                self.parent.setCentralWidget(central)
            else:
                layout = self.parent.centralWidget().layout()

            layout.addWidget(self.parent.graph_view)
