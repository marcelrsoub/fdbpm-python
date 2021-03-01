from PyQt5 import QtWidgets, QtGui

import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure

class MplWidget(FigureCanvasQTAgg):
    
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super(MplWidget, self).__init__(fig)

# --------------------------------------
class EasyPlotGUI(QtWidgets.QMainWindow):

    """Easy GUI with Matplotlib
    """
    
    def __init__(self, ui_python_file):
        self.app = QtWidgets.QApplication([])
        self.window_title="EasyPlotGUI"
    
        self.icon_path=None
        self.ui=ui_python_file.Ui_MainWindow()
        super().__init__()
        self.ui.setupUi(self) 


    def show_gui(self):
        """Shows the GUI.
        """

        try:
            self.ax=self.ui.MplWidget.axes
            self.draw=self.ui.MplWidget.draw
            self.update_graph()
        except AttributeError:
            pass

        self.update_interactivity()
        

        self.setWindowTitle(self.window_title)

        if self.icon_path!=None:
            
            self.setWindowIcon(QtGui.QIcon(self.icon_path))
        
        self.show()  # Show a window
        self.app.exec_()

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
        
        self.draw()