from PySide6.QtWidgets import QWidget, QVBoxLayout, QPushButton, QLabel, QFileDialog, QMessageBox
from PySide6.QtGui import QFont
from PySide6.QtCore import Qt, Signal

class FileSelector(QWidget):
    archivo_seleccionado = Signal(str)  # Señal que emite la ruta del archivo
    def __init__(self):
        super().__init__()

        self.layout = QVBoxLayout(self)
        self.layout.setAlignment(Qt.AlignCenter)  # Centrar todo el layout verticalmente y horizontalmente
        self.layout.setSpacing(15)  # Espacio entre widgets


        self.label = QLabel("Por favor, selecciona un archivo SDF para cargar.")
        font = QFont()
        font.setPointSize(16)
        font.setBold(True)     
        self.label.setFont(font)
        self.label.setAlignment(Qt.AlignCenter)
        self.layout.addWidget(self.label)

        self.open_button = QPushButton("Abrir archivo")
        self.open_button.setFont(font)
        self.layout.addWidget(self.open_button)

        self.selected_file = None

        self.open_button.clicked.connect(self.open_file_dialog)

    def open_file_dialog(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Seleccionar archivo SDF",
            "",
            "Archivos SDF (*.sdf);;Todos los archivos (*)"
        )
        if file_path:
            self.selected_file = file_path
            self.archivo_seleccionado.emit(file_path)  # Emitir la señal

