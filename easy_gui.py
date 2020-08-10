# -*- coding: utf-8 -*-
# @Author: Marcel Reis-Soubkovsky
# @Date:   2020-06-28 15:46:41
# @Last Modified by:   Marcel Reis-Soubkovsky
# @Last Modified time: 2020-07-08 17:37:59

from PySide2.QtWidgets import*
from PySide2.QtUiTools import QUiLoader
from PySide2.QtCore import Qt, QFile, QIODevice
from PySide2.QtWidgets import QApplication, QWidget
from PySide2.QtGui import QIcon
import sys

from matplotlib.backends.backend_qt5agg import (
        FigureCanvas, NavigationToolbar2QT as NavigationToolbar)

from matplotlib.figure import Figure
from matplotlib.widgets import Slider, Button, RadioButtons

import numpy as np
import random


class MplWidget(QWidget):
    
    def __init__(self, parent = None):
        
        QWidget.__init__(self, parent)
        
        self.canvas = FigureCanvas(Figure(constrained_layout=True, facecolor='#f0f0f0'))
        
        vertical_layout = QVBoxLayout()
        vertical_layout.addWidget(self.canvas)
        # vertical_layout.addWidget(NavigationToolbar(self.canvas, self))
        
        self.canvas.axes = self.canvas.figure.add_subplot(111)
        self.setLayout(vertical_layout)    

# ------------------ MainWidget ------------------
class EasyGUI(QWidget):

    """Easy GUI with matplotlib
    """
    
    def __init__(self):

        self.window_title="Python"
        self.ui_filepath="untitled.ui"
        self.icon_path=""

    def show_gui(self):
        app = QApplication(sys.argv)
        
        QWidget.__init__(self)

        loader = QUiLoader()
        loader.registerCustomWidget(MplWidget)
        ui_file = QFile(self.ui_filepath)
        if not ui_file.open(QIODevice.ReadOnly):
            print("Cannot open {}: {}".format(self.ui_filepath, ui_file.errorString()))
            sys.exit(-1)
        loader = QUiLoader()
        self.ui = loader.load(ui_file, None)
        ui_file.close()
        if not self.ui:
            print(loader.errorString())
            sys.exit(-1)

        

        

        grid_layout = QGridLayout()
        grid_layout.addWidget(self.ui)
        self.setLayout(grid_layout)

        self.setWindowTitle(self.window_title)
        # app.setWindowIcon(QIcon("beampy-logo.png"))
        

        self.update_interactivity()
        self.update_graph()

        self.setWindowIcon(QIcon("beampy-logo.png"))
        
        self.show()
        sys.exit(app.exec_())

    def update_interactivity(self):
        self.ui.pushButton_generate_random_signal.clicked.connect(self.update_graph)    # example, to be overwritten

    def update_graph(self):
        # example, to be overwritten
        fs = 500
        f = random.randint(1, 100)
        ts = 1/fs
        length_of_signal = 100
        t = np.linspace(0,1,length_of_signal)
        
        cosinus_signal = np.cos(2*np.pi*f*t)
        sinus_signal = np.sin(2*np.pi*f*t)

        self.ui.MplWidget.canvas.axes.clear()
        self.ui.MplWidget.canvas.axes.plot(t, cosinus_signal)
        self.ui.MplWidget.canvas.axes.plot(t, sinus_signal)
        self.ui.MplWidget.canvas.axes.legend(('cosinus', 'sinus'),loc='upper right')
        self.ui.MplWidget.canvas.axes.set_title('Cosinus - Sinus Signals')
        self.ui.MplWidget.canvas.draw()
        
if __name__ == "__main__":

    gui = EasyGUI()
    gui.window_title="TEST"
    gui.ui_filepath="untitled.ui"
    gui.show_gui()