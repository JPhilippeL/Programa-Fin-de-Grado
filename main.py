from PySide6.QtWidgets import QApplication
from ui.main_window import MainWindow
import sys

app = QApplication(sys.argv)
app.setStyle("Fusion")  # Establecer el estilo de la aplicación
app.setApplicationName("Editor Molecular")
window = MainWindow()
window.show()
sys.exit(app.exec())
