# -*- coding: utf-8 -*-
"""
Spyder Editor
This is a temporary script file.
https://stackoverflow.com/questions/52380572/handling-folds-in-spyder
"""

import sys
import qdarkstyle
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import ( QApplication, QComboBox, QDialog,
                                QDialogButtonBox, QFormLayout, QGridLayout, QGroupBox, QHBoxLayout,
                                QLabel, QLineEdit, QMenu, QMenuBar, QPushButton, QSpinBox, QTextEdit,
                                QVBoxLayout,QWidget )

from PyQt5.QtWidgets import QApplication, QWidget, QPushButton, QMessageBox
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import pyqtSlot

class Dialog(QDialog):
    NumGridRows = 3
    NumButtons = 4

    def __init__(self):
        super(Dialog, self).__init__()
        self.createFormGroupBox()

        buttonBox = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttonBox.accepted.connect(self.accept)
        buttonBox.rejected.connect(self.reject)

        mainLayout = QVBoxLayout()
        mainLayout.addWidget(self.formGroupBox)
        mainLayout.addWidget(buttonBox)
        self.setLayout(mainLayout)

        self.setWindowTitle("Form Layout - pythonspot.com")

    def createFormGroupBox(self):
        self.formGroupBox = QGroupBox("Form layout")
        layout = QFormLayout()
        layout.addRow(QLabel("Name:"), QLineEdit())
        layout.addRow(QLabel("Country:"), QComboBox())
        layout.addRow(QLabel("Age:"), QSpinBox())
        self.formGroupBox.setLayout(layout)

class MessageBox(QWidget):

    def __init__(self):
        super().__init__()
        self.title = 'PyQt5 messagebox - pythonspot.com'
        self.left = 10
        self.top = 10
        self.width = 320
        self.height = 200
        self.initUI()
        
    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        # self.buttonReply = QMessageBox.question(self, 'PyQt5 message', "Do you like PyQt5?", QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        self.buttonReply = QMessageBox.question(self, 'PyQt5 message', "Do you want to restart?", QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel, QMessageBox.Cancel)
        # self.show()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    mb = MessageBox()
    sys.exit(app.exec_())


class ACMGUI(object):
    def __init__(self):

        self.app = QApplication(sys.argv)
        self.mb = MessageBox()

        # # create the application and the main window
        # app = QApplication(sys.argv)
        # dialog = Dialog()
        #     # app = QtWidgets.QApplication(sys.argv)
        #     # window = QtWidgets.QMainWindow()

        # # setup stylesheet
        # app.setStyleSheet(qdarkstyle.load_stylesheet_pyqt5())
        # self.app = app

        # # run
        # # sys.exit(dialog.exec_())
        #     # window.show()
        #     # app.exec_()


if __name__ == '__main__':
    # create the application and the main window
    app = QApplication(sys.argv)
    dialog = Dialog()
    # app = QtWidgets.QApplication(sys.argv)
    # window = QtWidgets.QMainWindow()

    # setup stylesheet
    app.setStyleSheet(qdarkstyle.load_stylesheet_pyqt5())

    # run
    sys.exit(dialog.exec_())
    # window.show()
    # app.exec_()

