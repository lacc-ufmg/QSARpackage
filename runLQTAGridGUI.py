#!/usr/bin/env python

from PySide6.QtWidgets import QApplication, QMainWindow, QWidget, QVBoxLayout, QPushButton, QComboBox, QLineEdit, QFileDialog, QDoubleSpinBox, QLabel, QListWidget, QCheckBox, QHBoxLayout, QGroupBox, QFileDialog
from PySide6.QtCore import Qt
import sys
# Supondo que a lógica do seu programa esteja acessível como um módulo
from LQTAQSAR.LQTAGrid.grid_generate import *

class MainApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('Interface para LQTAGrid')
        self.setGeometry(100, 100, 800, 600)

        self.initUI()

    def initUI(self):
        centralWidget = QWidget(self)
        self.setCentralWidget(centralWidget)
        layout = QVBoxLayout(centralWidget)

        # Extensão
        self.extensionComboBox = QComboBox()
        self.extensionComboBox.addItems(['mol2', 'pdb', 'xyz'])
        layout.addWidget(QLabel('Extensão:'))
        layout.addWidget(self.extensionComboBox)

        # Coordenadas
        self.coordinatesLineEdits = [QLineEdit() for _ in range(3)]
        coordsLayout = QHBoxLayout()
        for le in self.coordinatesLineEdits:
            le.setPlaceholderText('Auto')
            coordsLayout.addWidget(le)
        coordsGroupBox = QGroupBox('Coordenadas:')
        coordsGroupBox.setLayout(coordsLayout)
        layout.addWidget(coordsGroupBox)

        # Dimensões
        self.dimensionsLineEdits = [QLineEdit() for _ in range(3)]
        dimsLayout = QHBoxLayout()
        for le in self.dimensionsLineEdits:
            le.setPlaceholderText('Auto')
            dimsLayout.addWidget(le)
        dimsGroupBox = QGroupBox('Dimensões:')
        dimsGroupBox.setLayout(dimsLayout)
        layout.addWidget(dimsGroupBox)

        # Átomo
        self.atomListWidget = QListWidget()
        self.atomListWidget.setSelectionMode(QListWidget.MultiSelection)
        self.atomListWidget.addItems(['NH3+', 'CH3COO-', 'Cl-', 'Na+'])
        layout.addWidget(QLabel('Átomo:'))
        layout.addWidget(self.atomListWidget)

        # Direório de Moléculas
        self.molsLineEdit = QLineEdit()
        self.molsButton = QPushButton('Selecionar Diretório')
        self.molsButton.clicked.connect(self.selectMolsDirectory)
        molsLayout = QHBoxLayout()
        molsLayout.addWidget(self.molsLineEdit)
        molsLayout.addWidget(self.molsButton)
        layout.addLayout(molsLayout)

        # Step
        self.stepSpinBox = QDoubleSpinBox()
        self.stepSpinBox.setValue(1.0)
        self.stepSpinBox.setMinimum(0.01)
        layout.addWidget(QLabel('Step:'))
        layout.addWidget(self.stepSpinBox)

        # Output
        self.outputLineEdit = QLineEdit()
        self.outputButton = QPushButton('Selecionar Arquivo de Saída')
        self.outputButton.clicked.connect(self.selectOutputFile)
        outputLayout = QHBoxLayout()
        outputLayout.addWidget(self.outputLineEdit)
        outputLayout.addWidget(self.outputButton)
        layout.addLayout(outputLayout)

        # Botão de execução
        self.runButton = QPushButton('Executar')
        self.runButton.clicked.connect(self.execute)
        layout.addWidget(self.runButton)

    def selectMolsDirectory(self):
        directory = QFileDialog.getExistingDirectory(self, "Selecionar Diretório de Moléculas")
        if directory:
            self.molsLineEdit.setText(directory)

    def selectOutputFile(self):
        outputFile, _ = QFileDialog.getSaveFileName(self, "Selecionar Arquivo de Saída", "", "CSV Files (*.csv)")
        if outputFile:
            self.outputLineEdit.setText(outputFile)

    def execute(self):
        extension = self.extensionComboBox.currentText()
        if any(not le.text() for le in self.coordinatesLineEdits):
            coordinates = None
        else:
            coordinates = [le.text() if le.text() else None for le in self.coordinatesLineEdits]
        if any(not le.text() for le in self.dimensionsLineEdits):
            dimensions = None
        else:
            dimensions = [le.text() if le.text() else None for le in self.dimensionsLineEdits]
        atom = [item.text() for item in self.atomListWidget.selectedItems()]
        mols = self.molsLineEdit.text()
        step = self.stepSpinBox.value()
        output = self.outputLineEdit.text()

        # Aqui você chamaria a função de cálculo do seu programa
        # Por exemplo: funcao_de_calculo(extension, coordinates, dimensions, atom, mols, step, output)
        print("Executando com:", extension, coordinates, dimensions, atom, mols, step, output)
        grid = GridGenerate(
        extension,
        coordinates,
        dimensions,
        atom,
        mols,
        step
        )
        grid.saveGrid(output)

        # Lembre-se de adaptar a chamada da função ao seu caso específico, incluindo a passagem de parâmetros corretos e o tratamento de valores opcionais.
        # Além disso, considere tratar exceções e fornecer feedback ao usuário sobre o sucesso da operação ou possíveis erros.

def main():
    app = QApplication(sys.argv)
    window = MainApp()
    window.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
t
