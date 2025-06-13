from PySide6.QtWidgets import QMenuBar, QFileDialog
from PySide6.QtGui import QAction
from core.sdf_parser import parse_sdf
from ui.graph_view import MoleculeGraphView
from core.sdf_parser import graph_to_mol
from rdkit import Chem
import networkx as nx

class MenuBar(QMenuBar):
    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent  # Referencia a MainWindow

        archivo_menu = self.addMenu("Archivo")

        nuevo_action = QAction("Nuevo", self)
        nuevo_action.triggered.connect(self.nuevo_archivo)
        archivo_menu.addAction(nuevo_action)


        cargar_action = QAction("Cargar", self)
        cargar_action.triggered.connect(self.cargar_archivo)
        archivo_menu.addAction(cargar_action)

        guardar_action = QAction("Guardar", self)
        guardar_action.triggered.connect(self.guardar_archivo)
        archivo_menu.addAction(guardar_action)  # sin funcionalidad por ahora

    def nuevo_archivo(self):
        self.parent.create_new_graph()

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

    def guardar_archivo(self):
        if not self.parent.graph_view:
            return  # No hay grafo cargado

        file_path, _ = QFileDialog.getSaveFileName(
            self.parent,
            "Guardar archivo SDF",
            "",
            "Archivos SDF (*.sdf)"
        )
        if not file_path:
            return

        # Convertir grafo a Mol y guardar
        mol = graph_to_mol(self.parent.graph_view.scene().graph)

        writer = Chem.SDWriter(file_path)
        writer.write(mol)
        writer.close()

