from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QFormLayout, QLineEdit,
    QPushButton, QFileDialog, QDialogButtonBox
)
import os

class ModelTestDialog(QDialog):
    last_model_path = ""
    last_sdf_path = ""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Seleccionar modelo y molécula")

        self.model_path_input = QLineEdit()
        self.model_path_input.setText(ModelTestDialog.last_model_path)
        self.model_browse_btn = QPushButton("Seleccionar...")
        self.model_browse_btn.clicked.connect(self.browse_model)

        self.sdf_path_input = QLineEdit()
        self.sdf_path_input.setText(ModelTestDialog.last_sdf_path)
        self.sdf_browse_btn = QPushButton("Seleccionar...")
        self.sdf_browse_btn.clicked.connect(self.browse_sdf)

        form_layout = QFormLayout()
        form_layout.addRow("Modelo (.pt):", self._with_button(self.model_path_input, self.model_browse_btn))
        form_layout.addRow("Molécula (.sdf):", self._with_button(self.sdf_path_input, self.sdf_browse_btn))

        layout = QVBoxLayout()
        layout.addLayout(form_layout)

        self.buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        layout.addWidget(self.buttons)

        self.setLayout(layout)

    def _with_button(self, line_edit, button):
        from PySide6.QtWidgets import QWidget, QHBoxLayout
        container = QWidget()
        hbox = QHBoxLayout(container)
        hbox.addWidget(line_edit)
        hbox.addWidget(button)
        hbox.setContentsMargins(0, 0, 0, 0)
        return container

    def browse_model(self):
        path, _ = QFileDialog.getOpenFileName(self, "Seleccionar modelo", "", "Modelos (*.pt)")
        if path:
            self.model_path_input.setText(path)

    def browse_sdf(self):
        path, _ = QFileDialog.getOpenFileName(self, "Seleccionar molécula", "", "Moléculas (*.sdf)")
        if path:
            self.sdf_path_input.setText(path)

    def get_paths(self):
        model_path = self.model_path_input.text()
        sdf_path = self.sdf_path_input.text()
        # Guardar para la proxima consulta
        ModelTestDialog.last_model_path = model_path
        ModelTestDialog.last_sdf_path = sdf_path
        return model_path, sdf_path

