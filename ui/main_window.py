from PySide6.QtWidgets import QMainWindow, QFileDialog, QWidget, QVBoxLayout, QMessageBox
from ui.graph_view import MoleculeGraphView
from core.sdf_parser import parse_sdf
from ui.file_selector import FileSelector
from ui.menu_bar import MenuBar
import networkx as nx

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Editor Molecular")
        self.resize(800, 600)

        # Contenedor central permanente
        self.central_widget = QWidget()
        self.central_layout = QVBoxLayout(self.central_widget)
        self.setCentralWidget(self.central_widget)

        self.menu_bar = MenuBar(self)
        self.setMenuBar(self.menu_bar)


        # Inicialmente, mostramos el selector de archivo
        self.file_selector = FileSelector()
        self.central_layout.addWidget(self.file_selector)

        # Conexión del botón "Abrir archivo"
        self.file_selector.open_button.clicked.connect(self.load_graph_if_selected)

        # Para luego reemplazarlo
        self.graph_view = None

    def create_new_graph(self):
        graph = nx.Graph()
        new_graph_view = MoleculeGraphView(graph)

        # Eliminar selector si existe
        if self.file_selector:
            self.central_layout.removeWidget(self.file_selector)
            self.file_selector.setParent(None)
            self.file_selector = None

        # Eliminar grafo anterior si existe
        if self.graph_view:
            self.central_layout.removeWidget(self.graph_view)
            self.graph_view.setParent(None)

        self.graph_view = new_graph_view
        self.central_layout.addWidget(self.graph_view)

        QMessageBox.information(
            self,
            "Consejo",
            "Haz clic derecho sobre el área para añadir un nodo."
        )


    def load_graph_if_selected(self):
        file_path = self.file_selector.selected_file
        if not file_path:
            return

        try:
            graph = parse_sdf(file_path)
        except Exception as e:
            QMessageBox.critical(self, "Error al cargar", f"No se pudo cargar el archivo:\n{str(e)}")
            return

        # Crear nueva vista de grafo
        new_graph_view = MoleculeGraphView(graph)

        # Reemplazar el selector de archivo por la vista del grafo
        self.central_layout.removeWidget(self.file_selector)
        self.file_selector.setParent(None)

        if self.graph_view:
            self.central_layout.removeWidget(self.graph_view)
            self.graph_view.setParent(None)

        self.graph_view = new_graph_view
        self.central_layout.addWidget(self.graph_view)
                
    def load_graph_from_file(self, file_path):
        try:
            graph = parse_sdf(file_path)
        except Exception as e:
            QMessageBox.critical(self, "Error al cargar", f"No se pudo cargar el archivo:\n{str(e)}")
            return

        new_graph_view = MoleculeGraphView(graph)

        # Si hay algo mostrado, lo quitamos
        if self.file_selector:
            self.central_layout.removeWidget(self.file_selector)
            self.file_selector.setParent(None)
            self.file_selector = None

        if self.graph_view:
            self.central_layout.removeWidget(self.graph_view)
            self.graph_view.setParent(None)

        self.graph_view = new_graph_view
        self.central_layout.addWidget(self.graph_view)
