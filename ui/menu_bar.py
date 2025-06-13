from PySide6.QtWidgets import QMenuBar, QFileDialog, QMessageBox
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
            self.parent.load_graph_from_file(file_path)


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

        # Convertir grafo a Mol
        try:
            mol = graph_to_mol(self.parent.graph_view.scene().graph)

            # Intentar sanitizar la molécula para detectar errores antes de guardar
            Chem.SanitizeMol(mol)  # Esto lanza excepciones si hay errores químicos

            # Intentar escribir el archivo
            writer = Chem.SDWriter(file_path)
            writer.write(mol)
            writer.close()

        except Exception as e:
            # Mostrar mensaje de error al usuario
            QMessageBox.critical(
                self.parent,
                "Error al guardar",
                f"No se pudo guardar la molécula:\n\n{str(e)}"
            )

