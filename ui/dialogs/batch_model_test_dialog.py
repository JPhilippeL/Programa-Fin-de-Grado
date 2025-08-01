from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QFormLayout, QLineEdit,
    QPushButton, QFileDialog, QDialogButtonBox, QWidget, QHBoxLayout
)

class BatchModelTestDialog(QDialog):
    last_model_path = ""
    last_sdf_dir = ""
    last_targets_file = ""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Testear modelo con conjunto de moléculas")

        # Inputs
        self.model_path_input = QLineEdit()
        self.model_path_input.setText(BatchModelTestDialog.last_model_path)
        self.model_browse_btn = QPushButton("Seleccionar...")
        self.model_browse_btn.clicked.connect(self.browse_model)

        self.sdf_dir_input = QLineEdit()
        self.sdf_dir_input.setText(BatchModelTestDialog.last_sdf_dir)
        self.sdf_browse_btn = QPushButton("Seleccionar...")
        self.sdf_browse_btn.clicked.connect(self.browse_sdf_dir)

        self.targets_file_input = QLineEdit()
        self.targets_file_input.setText(BatchModelTestDialog.last_targets_file)
        self.targets_browse_btn = QPushButton("Seleccionar...")
        self.targets_browse_btn.clicked.connect(self.browse_targets_file)

        # Form layout
        form_layout = QFormLayout()
        form_layout.addRow("Modelo (.pt):", self._with_button(self.model_path_input, self.model_browse_btn))
        form_layout.addRow("Directorio de SDFs:", self._with_button(self.sdf_dir_input, self.sdf_browse_btn))
        form_layout.addRow("Archivo de targets (.txt):", self._with_button(self.targets_file_input, self.targets_browse_btn))

        # Main layout
        layout = QVBoxLayout()
        layout.addLayout(form_layout)

        self.buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        layout.addWidget(self.buttons)

        self.setLayout(layout)

    def _with_button(self, line_edit, button):
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

    def browse_sdf_dir(self):
        path = QFileDialog.getExistingDirectory(self, "Seleccionar directorio de SDFs")
        if path:
            self.sdf_dir_input.setText(path)

    def browse_targets_file(self):
        path, _ = QFileDialog.getOpenFileName(self, "Seleccionar archivo de targets", "", "Texto (*.txt)")
        if path:
            self.targets_file_input.setText(path)

    def get_paths(self):
        model_path = self.model_path_input.text()
        sdf_dir = self.sdf_dir_input.text()
        targets_file = self.targets_file_input.text()

        # Guardar para la próxima vez
        BatchModelTestDialog.last_model_path = model_path
        BatchModelTestDialog.last_sdf_dir = sdf_dir
        BatchModelTestDialog.last_targets_file = targets_file

        return model_path, sdf_dir, targets_file
