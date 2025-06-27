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
        self.verificacion_al_guardar_activada = True  # Activada por defecto

        #Menu de Archivo    
        archivo_menu = self.addMenu("Archivo")

        nuevo_action = QAction("Nuevo", self)
        nuevo_action.triggered.connect(self.nuevo_archivo)
        archivo_menu.addAction(nuevo_action)

        cargar_action = QAction("Cargar", self)
        cargar_action.triggered.connect(self.cargar_archivo)
        archivo_menu.addAction(cargar_action)

        guardar_action = QAction("Guardar", self)
        guardar_action.triggered.connect(self.guardar_archivo)
        archivo_menu.addAction(guardar_action)

        # Menu de Verificación
        verificacion_menu = self.addMenu("Verificación")

        verificar_action = QAction("Verificar", self)
        verificar_action.triggered.connect(self.verificar_molecula)
        verificacion_menu.addAction(verificar_action)

        self.verificacion_guardar_action = QAction("Verificación al guardar", self)
        self.verificacion_guardar_action.setCheckable(True)
        self.verificacion_guardar_action.setChecked(True)
        self.verificacion_guardar_action.toggled.connect(self.toggle_verificacion_al_guardar)
        verificacion_menu.addAction(self.verificacion_guardar_action)

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

        sanitize = self.verificacion_al_guardar_activada
        mol = graph_to_mol(self.parent.graph_view.scene().graph, sanitize=sanitize)

        try:
            if sanitize:
                Chem.SanitizeMol(mol)

            writer = Chem.SDWriter(file_path)
            writer.write(mol)
            writer.close()
        except Exception as e:
            QMessageBox.critical(
                self.parent,
                "Error al guardar",
                f"No se pudo guardar la molécula:\n\n{str(e)}"
            )


    def verificar_molecula(self):
        if not self.parent.graph_view:
            QMessageBox.warning(
                self.parent,
                "Verificación",
                "No hay una molécula cargada para verificar."
            )
            return

        try:
            mol = graph_to_mol(self.parent.graph_view.scene().graph)
            Chem.SanitizeMol(mol)
            QMessageBox.information(
                self.parent,
                "Verificación Exitosa",
                "La molécula no contiene errores químicos detectables."
            )
        except Exception as e:
            QMessageBox.critical(
                self.parent,
                "Error de Verificación",
                f"Se detectaron errores químicos en la molécula:\n\n{str(e)}"
            )

    def toggle_verificacion_al_guardar(self, checked):
        self.verificacion_al_guardar_activada = checked

