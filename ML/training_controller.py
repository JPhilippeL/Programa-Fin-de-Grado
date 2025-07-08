# training_controller.py

from PySide6.QtCore import QObject, Signal, QThread
import logging
from ML.model_trainer import train_and_save_model

class TrainingController:
    def __init__(self, parent):
        self.parent = parent  # Puede ser MainWindow
        self.thread = None
        self.worker = None

    def entrenar(self, sdf_dir, target_file, modelo, epochs, batch_size, lr, valid_split, save_path):
        self.thread = QThread()
        self.worker = TrainerWorker(
            sdf_dir, target_file, modelo, epochs,
            save_path, batch_size=batch_size,
            lr=lr, valid_split=valid_split
        )

        self.worker.moveToThread(self.thread)

        # Conectar señales
        self.thread.started.connect(self.worker.run)
        self.worker.finished.connect(self.on_finished)
        self.worker.error.connect(self.on_error)
        self.worker.log.connect(self.parent.log)

        # Finalizar hilo correctamente
        self.worker.finished.connect(self.thread.quit)
        self.worker.finished.connect(self.worker.deleteLater)
        self.thread.finished.connect(self.thread.deleteLater)

        self.thread.start()

    def on_finished(self, ruta):
        from PySide6.QtWidgets import QMessageBox
        QMessageBox.information(self.parent, "Entrenamiento completo", f"Modelo guardado en: {ruta}")
        self.parent.log(f"Entrenamiento completo. Modelo en: {ruta}")

    def on_error(self, msg):
        from PySide6.QtWidgets import QMessageBox
        QMessageBox.critical(self.parent, "Error de entrenamiento", msg)
        self.parent.log(f"Error en entrenamiento: {msg}")

class TrainerWorker(QObject):
    log = Signal(str)
    finished = Signal(str)  # Ruta del modelo guardado
    error = Signal(str)

    def __init__(self, sdf_dir, target_file, modelo_nombre, epochs, save_path, batch_size=32, lr=0.001, valid_split=0.2):
        super().__init__()
        self.sdf_dir = sdf_dir
        self.target_file = target_file
        self.modelo_nombre = modelo_nombre
        self.epochs = epochs
        self.save_path = save_path
        self.batch_size = batch_size
        self.lr = lr
        self.valid_split = valid_split


    def run(self):
        try:
            path = train_and_save_model(
                sdf_dir=self.sdf_dir,
                target_file=self.target_file,
                modelo_nombre=self.modelo_nombre,
                epochs=self.epochs,
                save_path=self.save_path,
                batch_size=self.batch_size,
                lr=self.lr,
                valid_split=self.valid_split  # <- Nuevo
            )
            self.finished.emit(path)
        except Exception as e:
            logging.getLogger(__name__).exception("Error en entrenamiento")
            self.error.emit(str(e))

