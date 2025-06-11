from PySide6.QtWidgets import QMainWindow, QFileDialog, QWidget, QVBoxLayout
from ui.graph_view import MoleculeGraphView
from core.sdf_parser import parse_sdf

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Editor Molecular")
        self.resize(800, 600)

        # cargar sdf y grafo
        graph = parse_sdf("data/5RHD_ligand.sdf")

        # crear y mostrar grafo
        self.graph_view = MoleculeGraphView(graph)
        central = QWidget()
        layout = QVBoxLayout(central)
        layout.addWidget(self.graph_view)
        self.setCentralWidget(central)
