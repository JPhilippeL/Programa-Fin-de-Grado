# main.py
from PySide6.QtWidgets import QApplication
from PySide6.QtGui import QIcon
from ui.main_window import MainWindow
import sys

app = QApplication(sys.argv)
app.setStyle("Fusion")  # Establecer el estilo de la aplicación
app.setApplicationName("Editor Molecular")
app.setWindowIcon(QIcon("assets/icono.png"))
window = MainWindow()
window.show()
sys.exit(app.exec())
