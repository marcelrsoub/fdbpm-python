# -*- coding: utf-8 -*-
# @Author: Marcel Reis-Soubkovsky
# @Date:   2020-06-28 15:46:41
# @Last Modified by:   Marcel Reis-Soubkovsky
# @Last Modified time: 2020-07-08 18:04:03

# from PySide2.QtWidgets import *
from PySide2.QtUiTools import QUiLoader
from PySide2.QtCore import Qt, QFile, QIODevice
from PySide2.QtWidgets import QApplication, QWidget, QFileDialog, QVBoxLayout, QGridLayout, QMainWindow
from PySide2.QtGui import QIcon
import sys

from matplotlib.backends.backend_qt5agg import (
        FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
import matplotlib.pyplot as plt

from matplotlib.figure import Figure
from matplotlib.widgets import Slider, Button, RadioButtons
from cycler import cycler

import numpy as np


class MplWidget(QWidget):
    
    def __init__(self, parent = None):
        
        # super().__init__()
        QWidget.__init__(self, parent)
        
        self.mpl_canvas = FigureCanvas(Figure(constrained_layout=True, facecolor='#f0f0f0'))
        self.mpl_canvas.axes = self.mpl_canvas.figure.add_subplot(111)
        # self.mpl_canvas.setParent()
        
        vertical_layout = QVBoxLayout()
        vertical_layout.addWidget(self.mpl_canvas)
        self.setLayout(vertical_layout)
        # vertical_layout.addWidget(NavigationToolbar(self.mpl_canvas, self))       #uncomment in order to have the matplotlib tools
        

# --------------------------------------
class EasyGUI(QWidget):

    """Easy GUI with Matplotlib
    """
    
    def __init__(self):
        self.window_title="EasyPlotGUI"
        self.ui_filepath="test_example.ui"
        self.icon_path="logo.png"

        

    def show_gui(self):
        """This is the core of the execution.
        """
        app = QApplication([])
        
        QWidget.__init__(self)

        loader = QUiLoader()
        loader.registerCustomWidget(MplWidget)
        # self.ui = loader.load(designer_file, self)
        ui_file = QFile(self.ui_filepath)
        try:
            ui_file.open(QIODevice.ReadOnly)
        except:
            raise TypeError("Unable to open {}: {}".format(self.ui_filepath, ui_file.errorString()))
        # loader = QUiLoader()

        self.ui = loader.load(ui_file, None)
        ui_file.close()

        try:
            self.ax=self.ui.MplWidget.mpl_canvas.axes
            print("après")
            self.canvas=self.ui.MplWidget.mpl_canvas

            custom_cycler = (cycler(color=['#ed2124','#808081','c', 'm', 'y', 'k']))
            plt.rc('axes', prop_cycle=custom_cycler)
            self.update_graph()
        except AttributeError:
            # raise AttributeError("ERROR")
            pass

        if not self.ui:
            print(loader.errorString())
            sys.exit(-1)

        self.update_interactivity()
        

        self.setWindowTitle(self.window_title)

        

        grid_layout = QGridLayout()
        grid_layout.addWidget(self.ui)
        self.setLayout(grid_layout)

        self.QFileDialog=QFileDialog
        self.QMainWindow=QMainWindow

        if self.icon_path!=None:
            self.setWindowIcon(QIcon(self.icon_path))
        
        self.show()
        app.exec_()

    def update_interactivity(self):
        """This method adds the interactivity between the GUI elements and the code. An example of execution can be seen below. This method is supposed to be overwritten when EasyGUI is imported as a parent class for adding your own graphs.

        Check interactivity functions on PySide2 doc.

        """  
        # example, to be overwritten
        self.ui.pushButton_generate_random_signal.clicked.connect(self.update_graph)

    def update_graph(self):
        """This method generates and updates the graph. An example of execution can be seen below. This method is supposed to be overwritten when EasyGUI is imported as a parent class for adding your own graphs.

        """        
        # example, to be overwritten
        y = np.random.rand(100)
        x = np.linspace(0,1,100)

        self.ax.clear()
        self.ax.plot(x, y, label="Random")
        self.ax.legend()
        self.ax.set_title('Random Values')
        self.canvas.draw()
        
if __name__ == "__main__":
    # Test Example
    gui = EasyGUI()
    gui.window_title="Window Title"
    gui.ui_filepath="test_example.ui"
    gui.icon_path="./logo.png"
    gui.show_gui()

    # EXAMPLE OF USAGE
    # Slider changing the frequency of a sine wave.

    # class MyClass(EasyPlotGUI):
    #     def __init__(self):
    #         super().__init__()
    #         self.window_title="My GUI Name"
    #         self.ui_filepath="X:/xxxxx/xxxx/your_GUI.ui"
    #         self.icon_path="X:/xxxxx/xxxx/your_GUI_icon.png"

    #         #initialize Graph variables for first plot
    #         self.f=1
        
    #     def update_interactivity(self):
    #         self.ui.mySlider.valueChanged(self.change_frequency())
        
    #     def change_frequency(self):
    #         self.f=self.ui.mySlider.value()
    #         self.update_graph()
        
    #     def update_graph(self):
    #         x=np.linspace(0,1)
    #         y=np.sin(2*np.pi*self.f*x)

    #         self.ax.clear()
    #         self.ax.plot(x, y, label="Sine")
    #         self.ax.legend()
    #         self.ax.set_title('Sine Wave')
    #         self.canvas.draw()

    # #calling it
    # my_gui=MyClass()
    # my_gui.show_gui()


