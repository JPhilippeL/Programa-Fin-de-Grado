from PySide6.QtWidgets import QMenuBar, QFileDialog, QMessageBox, QInputDialog, QDialog
from PySide6.QtGui import QAction
from core.sdf_converter import graph_to_mol, save_graph_as_sdf
from ML.model_tester import cargar_y_predecir
import os
from rdkit import Chem
from ui.forms.train_config_dialog import TrainConfigDialog
import logging
logger = logging.getLogger(__name__)

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

        # Menu de IA
        ia_menu = self.addMenu("IA")

        # Entrenamiento de IA
        entrenar_action = QAction("Entrenar IA", self)
        entrenar_action.triggered.connect(self.entrenar_ia)
        ia_menu.addAction(entrenar_action)

        # Testeo de IA
        testeo_action = QAction("Testear IA", self)
        testeo_action.triggered.connect(self.testear_modelo)
        ia_menu.addAction(testeo_action)

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

        graph = self.parent.graph_view.scene().graph

        try:
            save_graph_as_sdf(graph, file_path)
            mensaje = f"Archivo guardado correctamente en: {file_path}"
            logger.info(mensaje)
        except Exception as e:
            QMessageBox.critical(
                self.parent,
                "Error al guardar",
                str(e)
            )


    def verificar_molecula(self):
        if not self.parent.graph_view:
            mensaje = "No hay una molécula cargada para verificar."
            QMessageBox.warning(self.parent, "Verificación", mensaje)
            return

        try:
            mol = graph_to_mol(self.parent.graph_view.scene().graph)
            Chem.SanitizeMol(mol)
            mensaje = "La molécula no contiene errores químicos detectables."
            logger.info(mensaje)
        except Exception as e:
            mensaje = f"Se detectaron errores químicos en la molécula:\n{str(e)}"
            #QMessageBox.critical(self.parent, "Error de Verificación", mensaje)
            logger.error(mensaje)

    def entrenar_ia(self):
        sdf_dir = QFileDialog.getExistingDirectory(self.parent, "Seleccionar carpeta con archivos SDF")
        if not sdf_dir:
            return

        target_file, _ = QFileDialog.getOpenFileName(
            self.parent, "Seleccionar archivo de propiedades", "", "Archivo de texto (*.txt)"
        )
        if not target_file:
            return

        # Mostrar formulario
        dialog = TrainConfigDialog(self)
        if dialog.exec_() != QDialog.Accepted:
            return

        config = dialog.get_values()
        if not config["save_name"]:
            QMessageBox.warning(self.parent, "Nombre inválido", "El nombre del archivo no puede estar vacío.")
            return

        save_path = f"modelos/{config['save_name']}.pt"

        self.parent.training_controller.entrenar(
            sdf_dir=sdf_dir,
            target_file=target_file,
            modelo=config["modelo"],
            epochs=config["epochs"],
            batch_size=config["batch_size"],
            lr=config["lr"],
            valid_split=config["valid_split"],
            save_path=save_path,
            hidden_dim=config["hidden_dim"],
            num_layers=config["num_layers"]
        )

    def testear_modelo(self):
        model_path, _ = QFileDialog.getOpenFileName(self.parent, "Seleccionar archivo de modelo (.pt)", "", "Modelos (*.pt)")
        if not model_path:
            return

        sdf_path, _ = QFileDialog.getOpenFileName(self.parent, "Seleccionar archivo SDF para predecir", "", "Archivos SDF (*.sdf)")
        if not sdf_path:
            return

        try:
            pred, target_name = cargar_y_predecir(model_path, sdf_path)

            model_name = os.path.basename(model_path)
            sdf_name = os.path.basename(sdf_path)

            msg = f"Predicción de '{target_name}' con el modelo '{model_name}' en la molécula '{sdf_name}': {pred:.4f}"
            QMessageBox.information(self.parent, "Predicción", msg)
            logger.info(msg)


        except Exception as e:
            QMessageBox.critical(self.parent, "Error en predicción", f"No se pudo realizar la predicción:\n\n{str(e)}")
            logger.error(f"Error en testear modelo: {str(e)}")




