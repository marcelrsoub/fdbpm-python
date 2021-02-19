import numpy as np
import time as t

from PySide2.QtWidgets import *

from cycler import cycler

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from functools import partial
import copy

from easy_gui import EasyGUI

from numba import jit


class FdBpm(EasyGUI):
    def __init__(self):
        super().__init__()
        self.slider_realtime = True
        self.window_title = "FD-BPM"
        self.ui_filepath = r"interface.ui"
        self.icon_path = r"fdbpm-logo.png"

    def create_space(self):
        # window properties
        if not hasattr(self, 'NUM_SAMPLES') or not hasattr(self, 'LENGTH'):
            self.NUM_SAMPLES = int(100)
            self.LENGTH = int(1E3)   # set to 1E3 for better precision
        self.cmap = "jet"
        # Physical Properties
        if not hasattr(self, 'L'):
            self.L = 240E-6
        if not hasattr(self, 'dy'):
            self.dy = 1e-6  # set to 1E6
        if not hasattr(self, 'l_ambda'):
            self.l_ambda = 0.8e-6

        self.dx = self.L/self.NUM_SAMPLES
        self.x = np.linspace(-self.L/2, +self.L/2, self.NUM_SAMPLES)

        self.light = np.zeros(self.x.shape)

    def create_source(self, waist=8e-6, offset=0, plotOn=True):

        self.light = self.gauss(self.x-offset, waist)

        self.light_offset = offset

        if plotOn:
            plt.plot(self.x, self.light)
            plt.show()

    def gauss_light(self, fwhm=20E-6, offset=0):
        """
        Create a gaussian beam in amplitude.

        :math:`E = e^{-((x-x_0)/w)^{2P}}`

        The waist is defined as fwhm/sqrt(2*log(2)) and correspond to the 1/e
        field value and 1/:math:`e^2` intensity value.

        Parameters
        ----------
        fwhm : float
            Full width at half maximum (for intensity not amplitude) (µm).
        offset : float, optional
            Light offset from center in µm. 0 by default.

        Returns
        -------
        field : array
            Amplitude values over x in µm.

        Notes
        -----
        This methods uses the x and dist_x variables defined in :class:`Bpm`.
        """
        if not hasattr(self, 'fwhm'):
            self.fwhm = fwhm
        # such as I=1/e^2 in intensity
        spot_size = self.fwhm / np.sqrt(2 * np.log(2))
        if spot_size != 0:
            field = np.exp(-(self.x / spot_size)**2)
            field = np.roll(field, int(round(offset / self.dx)))
        else:
            field = 0 * self.x  # Avoid division by zero error
        self.light_offset = offset
        self.light = field
        return field

    def gauss(self, x, w0):
        return np.array(np.exp(-(x/w0)**2))

    def create_guides(self, width=25E-6, offset=0, plotOn=False):
        if not hasattr(self, 'n_env'):
            self.n_env = 1.004
        if not hasattr(self, 'guides'):
            self.guides = np.ones((self.NUM_SAMPLES,))*self.n_env
        if not hasattr(self, 'dn'):
            self.dn = 0.058

        mask = np.logical_and(self.x > -width/2+offset,
                              self.x < width/2+offset)
        self.guides[mask] += self.dn  # square guide
        # avg_guide = (np.max(guides)+np.min(guides))/2
        # avg_guide = (np.min(guides))
        self.avg_guide = (np.mean(self.guides))

        if plotOn:
            plt.plot(self.x, self.guides)
            plt.show()

        return (self.guides, self.avg_guide)

    def make_tri_matrix(self):

        if not hasattr(self, 'guides'):
            self.n_env = 1.004
            self.guides = np.ones((self.NUM_SAMPLES,))*self.n_env
            self.avg_guide = self.guides

        k0 = 2*np.pi/self.l_ambda
        # k = k0*self.guides
        k = k0/self.guides
        k_bar = k*self.avg_guide
        # k_bar = k/self.avg_guide

        self.h = self.dy
        self.ro = self.dy/(self.dx**2)
        self.A = 1j/(2*k_bar)
        self.B = 1j*(k**2-k_bar**2)/(2*k_bar)
        a = -self.ro*self.A
        b = 2*(1+self.ro*self.A)-self.h*self.B

        tridiag_matrix = np.zeros(
            (self.NUM_SAMPLES, self.NUM_SAMPLES), dtype='complex_')

        index = np.arange(self.NUM_SAMPLES)

        index0 = index[:-1]
        tridiag_matrix[index0, index0+1] = a[index0]

        index1 = index[1:]
        tridiag_matrix[index1, index1-1] = a[index1]

        tridiag_matrix[index, index] = b[index]

        # TODO: scipy class sparse matrices

        self.tridiag_matrix = tridiag_matrix

    @staticmethod
    @jit(forceobj=True)
    def core_propag_jit(A,B, ro, light, tridiag_matrix, h,LENGTH,NUM_SAMPLES):

        propag = np.zeros((int(LENGTH), int(NUM_SAMPLES)))
        d = np.zeros((NUM_SAMPLES), dtype='complex_')

        index = np.arange(1, NUM_SAMPLES-1)

        # cannot be transferred to numpy cause of np.linalg.solve in between lines
        for n in range(int(LENGTH)):
            d[index] = (2*(1-ro*A[index])+h*B[index])*light[index] + \
                ro*A[index]*(light[index-1]+light[index+1])
            d[0] = (2*(1-ro*A[0])+h*B[0])*light[0]+ro*A[0]*(light[1])
            d[-1] = (2*(1-ro*A[-1])+h*B[-1])*light[-1]+ro*A[-1]*(light[-2])

            # light = np.linalg.lstsq(tridiag_matrix,d)[0]
            light = np.linalg.solve(tridiag_matrix, d)  # fastest option
            # propag[n,:] = np.abs(light)
            propag[n, :] = (light*light.conjugate()).real
            # propag[n,:] = np.real(light) # Oscillations
            # propag[n,:] = np.angle(light) # interesting
        return propag

    def calculate_propagation(self, plotOn=True, timing=True):

        if (not hasattr(self, 'tridiag_matrix')):
            self.make_tri_matrix()

        if(timing == True):
            start = t.time()

        propag = self.core_propag_jit(self.A,self.B, self.ro, self.light, self.tridiag_matrix, self.h,self.LENGTH,self.NUM_SAMPLES)

        if(timing == True):
            end = t.time()
            print("Time elapsed:%0.4f seconds" % (end-start))

        self.propag = propag

        if(plotOn == True):
            propag_img = propag[::int(
                np.floor(self.LENGTH/self.NUM_SAMPLES)), :]

            fig, ax = plt.subplots()
            ax.imshow(propag_img, cmap=self.cmap, interpolation='bilinear', extent=[
                      -self.L/2*1E6, +self.L/2*1E6, self.LENGTH*self.dy*1E6, 0], aspect='auto')
            ax.set_xlabel(r"x ($\mu$m)")
            ax.set_ylabel(r"Length ($\mu$m)")

            plt.show()

        return self.propag

    def plot_moving_source(self):
        plt.close('all')

        fig, ax = plt.subplots()
        plt.subplots_adjust(left=0.25, bottom=0.25)
        if int(np.floor(self.LENGTH/self.NUM_SAMPLES)) > 0:
            propag_img = self.calculate_propagation(
                plotOn=False)[::int(np.floor(self.LENGTH/self.NUM_SAMPLES)), :]
        else:
            propag_img = self.calculate_propagation(plotOn=False)
        img = ax.imshow(propag_img, cmap=self.cmap, interpolation='bilinear', extent=[
                        -self.L/2*1E6, +self.L/2*1E6, self.LENGTH*self.dy*1E3, 0], aspect='auto')
        ax.set_xlabel(r"x ($\mu$m)")
        ax.set_ylabel(r"Length (mm)")
        ax.margins(x=0)

        axcolor = 'lightgoldenrodyellow'
        ax_position = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
        plt.plot(self.x, self.guides, 'orange')

        slider_position = Slider(ax_position, 'Light Offset', -self.L/2, +self.L/2,
                                 valinit=self.light_offset, valstep=self.dx*1E-2, valfmt="%2.0e")

        def update(val):
            pos = np.float(slider_position.val)
            # self.create_source(offset=pos,plotOn=False)
            self.gauss_light(offset=pos)
            img.set_data(self.calculate_propagation(plotOn=False))
            fig.canvas.draw_idle()

        slider_position.on_changed(update)

        plt.show()

    def reset_variables(self):
        # pass
        delattr(self, 'NUM_SAMPLES')
        delattr(self, 'LENGTH')
        delattr(self, 'dy')
        delattr(self, 'dx')
        delattr(self, 'l_ambda')
        delattr(self, 'L')
        delattr(self, 'fwhm')
        delattr(self, 'dn')
        delattr(self, 'light')
        delattr(self, 'light_offset')
        delattr(self, 'n_env')
        delattr(self, 'x')
        delattr(self, 'guides')
        delattr(self, 'avg_guide')
        delattr(self, 'h')
        delattr(self, 'ro')
        delattr(self, 'A')
        delattr(self, 'B')
        delattr(self, 'tridiag_matrix')

    def update_interactivity(self):
        self.ui.doubleSpinBox_lambda.setValue(self.l_ambda*1E6)
        self.ui.spinBox_propLength.setValue(self.LENGTH*self.dy*1E6)
        self.ui.spinBox_windowSize.setValue(self.L*1E6)
        self.ui.spinBox_numSamples.setValue(self.NUM_SAMPLES)

        self.ui.SliderPos.setMinimum(-self.L/2*1E6)
        self.ui.SliderPos.setMaximum(+self.L/2*1E6)
        self.ui.SliderPos.setValue(self.light_offset*1E6)
        self.ui.label_slider.setText(
            str(np.format_float_positional(self.light_offset*1E6, 3))+" um")

        self.ui.actionFree_Space_Propagation.triggered.connect(
            partial(self.examples, 'free_space', chart_update=True))
        self.ui.actionOne_Waveguide.triggered.connect(
            partial(self.examples, 'one_guide', chart_update=True))
        self.ui.actionThree_Waveguides.triggered.connect(
            partial(self.examples, '3_guides', chart_update=True))

        # self.ui.SliderPos.setTickInterval(self.dx)

        def update_label():
            pos = self.ui.SliderPos.value()
            self.ui.label_slider.setText(str(pos)+" um")

        def update_input_values():
            dy = copy.deepcopy(self.dy)
            fwhm = copy.deepcopy(self.fwhm)
            self.reset_variables()
            self.dy = dy
            self.fwhm = fwhm
            self.l_ambda = 1E-6*np.float(self.ui.doubleSpinBox_lambda.value())
            self.LENGTH = 1/self.dy * \
                np.float(self.ui.spinBox_propLength.value()*1E-6)
            self.L = 1E-6*np.float(self.ui.spinBox_windowSize.value())
            self.NUM_SAMPLES = int(self.ui.spinBox_numSamples.value())
            self.create_space()
            self.gauss_light()
            self.create_guides()
            self.ui.SliderPos.setValue(self.light_offset*1E6)
            self.ui.SliderPos.setMinimum(-self.L/2*1E6)
            self.ui.SliderPos.setMaximum(+self.L/2*1E6)
            self.update_graph()
            # FIXME: update doesn't refresh guides

        def changed_value():
            pos = self.ui.SliderPos.value()
            self.gauss_light(offset=pos*1E-6)
            self.update_graph(slider=True)
            self.ui.label_slider.setText(str(pos)+" um")

        self.ui.pushButton.clicked.connect(update_input_values)

        if self.slider_realtime == True:
            self.ui.SliderPos.valueChanged.connect(changed_value)

        else:
            self.ui.SliderPos.sliderReleased.connect(changed_value)
            self.ui.SliderPos.valueChanged.connect(update_label)

    def examples(self, name, chart_update=False):
        if chart_update:
            self.reset_variables()
        if name == 'free_space':

            self.NUM_SAMPLES = 101
            self.LENGTH = 1E2
            self.dy = 1E-4
            self.l_ambda = 1.5E-6
            self.L = 1000E-6

            self.create_space()
            self.gauss_light(fwhm=20E-6)
            self.create_guides(width=0)

        elif name == 'one_guide':
            self.NUM_SAMPLES = 101
            self.LENGTH = 1E2
            self.dy = 1E-4
            self.l_ambda = 1.5E-6
            self.L = 500E-6

            self.create_space()
            self.gauss_light()
            self.dn = 0.058
            self.n_env = 1.004
            self.create_guides()
        elif name == '3_guides':
            self.NUM_SAMPLES = 101
            self.LENGTH = 500
            self.dy = 1E-5
            self.l_ambda = 1.55E-6
            self.L = 60E-6
            self.create_space()
            plotOn = True
            offset = 14E-6
            self.gauss_light(fwhm=12E-6, offset=-offset+1E-7)
            self.create_guides(width=10E-6, offset=-offset)
            self.create_guides(width=10E-6)
            self.create_guides(width=10E-6, offset=offset, plotOn=False)

        if chart_update:
            self.ui.doubleSpinBox_lambda.setValue(self.l_ambda*1E6)
            self.ui.spinBox_propLength.setValue(self.LENGTH*self.dy*1E6)
            self.ui.spinBox_windowSize.setValue(self.L*1E6)
            self.ui.spinBox_numSamples.setValue(self.NUM_SAMPLES)
            self.ui.SliderPos.setMinimum(-self.L/2*1E6)
            self.ui.SliderPos.setMaximum(+self.L/2*1E6)
            self.ui.SliderPos.setValue(self.light_offset*1E6)
            self.update_graph(slider=False)

    def update_graph(self, slider=False):

        canvas = self.ui.MplWidget.mpl_canvas
        ax = self.ui.MplWidget.mpl_canvas.axes  # no canvas in the object

        if int(np.floor(self.LENGTH/self.NUM_SAMPLES)) > 0:
            propag_img = self.calculate_propagation(
                plotOn=False)[::int(np.floor(self.LENGTH/self.NUM_SAMPLES)), :]
        else:
            propag_img = self.calculate_propagation(plotOn=False)

        if slider == False:
            ax.clear()
            self.img = ax.imshow(propag_img, cmap=self.cmap, interpolation='bilinear', extent=[
                                 -self.L/2*1E6, +self.L/2*1E6, self.LENGTH*self.dy*1E3, 0], aspect='auto')
            ax.set_xlabel(r"x ($\mu$m)")
            ax.set_ylabel(r"Length (mm)")
        else:
            # ax.clear()
            self.img.set_data(propag_img)
        plt.tight_layout(pad=0.)
        canvas.draw()


if __name__ == "__main__":
    # %matplotlib qt
    fd = FdBpm()

    fd.examples('one_guide')
    fd.show_gui()
