from PySide6.QtWidgets import QMainWindow, QFileDialog, QWidget, QVBoxLayout
from ui.graph_view import MoleculeGraphView
from core.sdf_parser import parse_sdf
from ui.file_selector import FileSelector
from PySide6.QtWidgets import QMessageBox

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Editor Molecular")
        self.resize(800, 600)

        self.file_selector = FileSelector()
        self.setCentralWidget(self.file_selector)

        self.file_selector.open_button.clicked.connect(self.load_graph_if_selected)

    def load_graph_if_selected(self):
        file_path = self.file_selector.selected_file
        if not file_path:
            return

        try:
            graph = parse_sdf(file_path)
        except Exception as e:
            QMessageBox.critical(self, "Error al cargar", f"No se pudo cargar el archivo:\n{str(e)}")
            return

        # Reemplazar el widget selector por la vista del grafo
        self.graph_view = MoleculeGraphView(graph)
        self.setCentralWidget(self.graph_view)
