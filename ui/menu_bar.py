from PySide6.QtWidgets import QMenuBar, QFileDialog, QMessageBox, QInputDialog
from PySide6.QtGui import QAction
from core.sdf_parser import parse_sdf
from ui.graph_view import MoleculeGraphView
from core.sdf_parser import graph_to_mol
from ML.torch_data_processing import read_targets, load_data_from_sdf, create_dataloader
from ML.torch_model_trainer import create_model, train
import torch
from rdkit import Chem
import os

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

        mol = graph_to_mol(self.parent.graph_view.scene().graph)

        try:
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

    def entrenar_ia(self):
        sdf_dir = QFileDialog.getExistingDirectory(self.parent, "Seleccionar carpeta con archivos SDF")
        if not sdf_dir:
            return

        target_file, _ = QFileDialog.getOpenFileName(self.parent, "Seleccionar archivo de propiedades", "", "Archivo de texto (*.txt)")
        if not target_file:
            return

        modelos = ["GIN", "GINE", "GAT"]
        modelo, ok = QInputDialog.getItem(self.parent, "Seleccionar modelo", "Modelo:", modelos, 0, False)
        if not ok:
            return

        epochs, ok = QInputDialog.getInt(self.parent, "Épocas de entrenamiento", "Número de épocas:", 20, 1, 10000)
        if not ok:
            return

        save_name, ok = QInputDialog.getText(self.parent, "Guardar modelo", "Nombre del archivo de modelo:")
        if not ok or not save_name:
            return
        save_path = f"modelos/{save_name}.pt"

        try:
            # Carga y preprocesa datos
            target_dict = read_targets(target_file)
            data_list = load_data_from_sdf(sdf_dir, target_dict)
            dataloader = create_dataloader(data_list, batch_size=32)

            device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

            # Crea el modelo elegido
            input_dim = 1      # Por tu mol_to_graph_data_obj: 1 feature por nodo (atomic number)
            edge_dim = 1       # bond type float
            model = create_model(modelo, input_dim, edge_dim)

            # Entrena el modelo
            train(model, dataloader, device, epochs=epochs, lr=0.001)

            checkpoint = {
                'model_state_dict': model.state_dict(),
                'model_type': modelo,
                'input_dim': input_dim,
                'edge_dim': edge_dim,
                'epochs_trained': epochs,
                # agrega otros parámetros que uses
            }

            os.makedirs(os.path.dirname(save_path), exist_ok=True)
            torch.save(checkpoint, save_path)

            QMessageBox.information(self.parent, "Entrenamiento completo", f"Modelo guardado en:\n{save_path}")

        except Exception as e:
            QMessageBox.critical(self.parent, "Error durante el entrenamiento", f"Ha ocurrido un error:\n\n{str(e)}")

    def testear_modelo(self):
        # Seleccionar archivo modelo .pt
        checkpoint_path, _ = QFileDialog.getOpenFileName(self.parent, "Seleccionar archivo de modelo (.pt)", "", "Modelos (*.pt)")
        if not checkpoint_path:
            return

        # Seleccionar archivo SDF
        sdf_path, _ = QFileDialog.getOpenFileName(self.parent, "Seleccionar archivo SDF para predecir", "", "Archivos SDF (*.sdf)")
        if not sdf_path:
            return

        try:
            from ML.model_tester import cargar_y_predecir
            pred = cargar_y_predecir(checkpoint_path, sdf_path)
            QMessageBox.information(self.parent, "Predicción", f"El target predicho es: {pred:.4f}")
            # en consola tmb
            print(f"Predicción para {sdf_path}: {pred:.4f}")
        except Exception as e:
            QMessageBox.critical(self.parent, "Error en predicción", f"No se pudo realizar la predicción:\n{str(e)}")
            # Q salga el error en la consola
            print(f"Error en testear modelo: {str(e)}")



