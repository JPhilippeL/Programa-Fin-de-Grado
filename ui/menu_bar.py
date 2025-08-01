from PySide6.QtWidgets import QMenuBar, QFileDialog, QMessageBox, QInputDialog, QDialog
from PySide6.QtGui import QAction
from core.sdf_converter import graph_to_mol, save_graph_as_sdf
from ML.model_tester import cargar_y_predecir
import os
from rdkit import Chem
from ui.dialogs.train_config_dialog import TrainConfigDialog
from ui.dialogs.model_test_dialog import ModelTestDialog
from ui.dialogs.batch_model_test_dialog import BatchModelTestDialog
from ML.model_tester import test_model_on_directory
import traceback
import torch
import logging
from ui.dialogs.image_dialog import ImageDialog
logger = logging.getLogger(__name__)

class MenuBar(QMenuBar):
    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent  # Referencia a MainWindow

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

        # Testeo de IA con un solo SDF
        testeo_action = QAction("Predecir SDF", self)
        testeo_action.triggered.connect(self.testear_modelo)
        ia_menu.addAction(testeo_action)

        # Testeo de IA con múltiples SDF
        testeo_batch_action = QAction("Testear IA", self)
        testeo_batch_action.triggered.connect(self.testear_modelo_en_batch)
        ia_menu.addAction(testeo_batch_action)

        # Consultar parámetros modelo
        consultar_params_action = QAction("Consultar modelo", self)
        consultar_params_action.triggered.connect(self.consultar_parametros_modelo)
        ia_menu.addAction(consultar_params_action)

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
            return

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
            self.parent.load_graph_from_file(file_path)  # Recargar el grafo guardado
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
            logger.error(mensaje)

    def entrenar_ia(self):
        dialog = TrainConfigDialog(self)
        if dialog.exec_() != QDialog.Accepted:
            return

        config = dialog.get_values()

        # Validaciones básicas
        if not config["sdf_dir"] or not os.path.isdir(config["sdf_dir"]):
            QMessageBox.warning(self.parent, "Error", "Debes seleccionar un directorio válido con archivos SDF.")
            return

        if not config["target_file"] or not os.path.isfile(config["target_file"]):
            QMessageBox.warning(self.parent, "Error", "Debes seleccionar un archivo .txt válido con los targets.")
            return

        if not config["save_name"]:
            QMessageBox.warning(self.parent, "Nombre inválido", "El nombre del archivo no puede estar vacío.")
            return

        save_path = f"modelos/{config['save_name']}.pt"

        # Ejecutar entrenamiento
        self.parent.training_controller.entrenar(
            sdf_dir=config["sdf_dir"],
            target_file=config["target_file"],
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
        dialog = ModelTestDialog(self.parent)
        if dialog.exec():
            model_path, sdf_path = dialog.get_paths()
            try:
                pred, target_name = cargar_y_predecir(model_path, sdf_path)

                model_name = os.path.basename(model_path)
                sdf_name = os.path.basename(sdf_path)
                msg = f"Predicción de '{target_name}' con el modelo '{model_name}' en la molécula '{sdf_name}': {pred:.4f}"
                logger.info(msg)

            except Exception as e:
                QMessageBox.critical(self.parent, "Error en predicción", f"No se pudo realizar la predicción:\n\n{str(e)}")
                logger.error(f"Error en testear modelo: {str(e)}")

    def testear_modelo_en_batch(self):

        dialog = BatchModelTestDialog(self.parent)
        if dialog.exec():
            model_path, sdf_dir, targets_file = dialog.get_paths()

            try:
                # Definir nombre del archivo de salida de predicciones
                model_name = os.path.splitext(os.path.basename(model_path))[0]

                # Ejecutar función de testeo
                test_model_on_directory(model_path, sdf_dir, targets_file)

                # Mostrar scatter plot
                plot_path = os.path.join("Resultados", f"scatter_plot_{model_name}.png")
                self.image_dialog = ImageDialog(plot_path, self.parent)
                self.image_dialog.show()



            except Exception as e:
                logger.error("Error en testeo por lotes:\n" + traceback.format_exc())
                QMessageBox.critical(self.parent, "Error en testeo", f"No se pudo ejecutar el testeo:\n\n{str(e)}")


    def consultar_parametros_modelo(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self.parent,
            "Seleccionar archivo de modelo (.pt)",
            "",
            "Modelos (*.pt)"
        )
        if not file_path:
            return

        try:
            checkpoint = torch.load(file_path, map_location='cpu')
            info = (
                f"Modelo: {checkpoint.get('model_type', 'Desconocido')}\n"
                f"\t\tÉpocas entrenadas: {checkpoint.get('epochs_trained', 'Desconocido')}\n"
                f"\t\tTarget: {checkpoint.get('target_name', 'Desconocido')}\n"
                f"\t\tHidden dim: {checkpoint.get('hidden_dim', 'Desconocido')}\n"
                f"\t\tNúmero de capas: {checkpoint.get('num_layers', 'Desconocido')}\n"
                f"\t\tBatch size: {checkpoint.get('batch_size', 'Desconocido')}\n"
                f"\t\tLearning rate: {checkpoint.get('learning_rate', 'Desconocido')}\n"
                f"\t\tValid split: {checkpoint.get('valid_split', 'Desconocido')}"
            )

            logger.info(info)

        except Exception as e:
            logger.error(f"Error al consultar parámetros del modelo: {str(e)}")






