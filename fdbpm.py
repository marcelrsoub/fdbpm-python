import numpy as np
import time as t
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

class FdBpm:
    def __init__(self):
        # window properties
        self.NUM_SAMPLES = np.int(500)
        self.LENGTH = np.int(1E3)   # set to 1E3 for better precision

    def create_space(self):
        # Physical Properties
        if not hasattr(self,'L'):
            self.L=240E-6
        if not hasattr(self,'dy'):
            self.dy = 1e-6  #set to 1E6
        if not hasattr(self,'l_ambda'):
            self.l_ambda = 0.8e-6

        self.dx = self.L/self.NUM_SAMPLES
        self.x = np.linspace(-self.L/2, +self.L/2, self.NUM_SAMPLES)
        
        self.light=np.zeros(self.x.shape)

    def create_source(self,waist=8e-6,offset=0,plotOn=True):
        
        self.light=self.gauss(self.x-offset,waist)

        self.light_offset=offset

        if plotOn:
            plt.plot(self.x,self.light)
            plt.show()
        

    def gauss(self,x,w0):
        return np.array(np.exp(-(x/w0)**2))

    def create_guides(self,width=8E-6,offset=0,plotOn=False):
        n_env=1.0405
        if not hasattr(self,'guides'):
            self.guides = np.ones((self.NUM_SAMPLES,))*n_env
        
        dn=0.03
        mask=np.logical_and(self.x>-width/2+offset ,self.x<width/2+offset)
        self.guides[mask]+=dn  #square guide
        # avg_guide = (np.max(guides)+np.min(guides))/2
        # avg_guide = (np.min(guides))
        self.avg_guide = (np.mean(self.guides))

        
        if plotOn:
            plt.plot(self.x,self.guides)
            plt.show()

        return (self.guides,self.avg_guide)

    def make_tri_matrix(self):

        k0 = 2*np.pi/self.l_ambda
        k = k0/self.guides
        k_bar = k*self.avg_guide

        self.h = self.dy
        self.ro = self.dy/(self.dx**2)
        self.A = 1j/(2*k_bar)
        self.B = 1j*(k**2-k_bar**2)/(2*k_bar)
        a = -self.ro*self.A
        b = 2*(1+self.ro*self.A)-self.h*self.B

        tridiag_matrix = np.zeros((self.NUM_SAMPLES,self.NUM_SAMPLES),dtype='complex_')

        index=np.arange(self.NUM_SAMPLES)

        index0=index[:-1]
        tridiag_matrix[index0,index0+1]=a[index0]

        index1=index[1:]
        tridiag_matrix[index1,index1-1]=a[index1]

        tridiag_matrix[index,index]=b[index]

        self.tridiag_matrix=tridiag_matrix

    def calculate_propagation(self, plotOn=True, timing=True):

        if not hasattr(self,'tridiag_matrix'):
            self.make_tri_matrix()

        if timing:
            start=t.time()

        (A,B,ro,light,tridiag_matrix,h)=(self.A,self.B,self.ro,self.light,self.tridiag_matrix,self.h)

        propag = np.zeros((np.int(self.LENGTH),np.int(self.NUM_SAMPLES)))
        d = np.zeros((self.NUM_SAMPLES),dtype='complex_')

        index=np.arange(1,self.NUM_SAMPLES-1)

        for n in range(np.int(self.LENGTH)):         # cannot be transferred to numpy cause of np.linalg.solve in between lines
            d[index]=(2*(1-ro*A[index])+h*B[index])*light[index]+ro*A[index]*(light[index-1]+light[index+1])
            d[0] = (2*(1-ro*A[0])+h*B[0])*light[0]+ro*A[0]*(light[1])
            d[-1] = (2*(1-ro*A[-1])+h*B[-1])*light[-1]+ro*A[-1]*(light[-2])

            # self.modo = np.linalg.lstsq(tridiag_matrix,d)[0]
            light = np.linalg.solve(tridiag_matrix,d) #fastest option
            propag[n,:] = np.abs(light)
            # propag[n,:] = light.real  # E Field

        if timing:
            end = t.time()
            print("Time elapsed:%0.4f seconds"%(end-start))

        self.propag=propag

        if plotOn:
            propag_img=propag[::np.int(np.floor(self.LENGTH/self.NUM_SAMPLES)),:]
            
            fig,ax=plt.subplots()
            ax.imshow(propag_img,cmap="plasma",interpolation='bilinear',extent=[-self.L/2*1E6,+self.L/2*1E6,self.LENGTH*self.dy*1E6,0],aspect='auto')
            ax.set_xlabel(r"x ($\mu$m)")
            ax.set_ylabel(r"Length ($\mu$m)")

            plt.show()

        return self.propag

    def plot_moving_source(self):
        plt.close('all')

        fig, ax = plt.subplots()
        plt.subplots_adjust(left=0.25, bottom=0.25)
        img=ax.imshow(self.calculate_propagation(plotOn=False),cmap="jet",interpolation='bilinear',extent=[-self.L/2*1E6,+self.L/2*1E6,self.LENGTH*self.dy*1E6,0],aspect='auto')
        ax.set_xlabel(r"x ($\mu$m)")
        ax.set_ylabel(r"Length ($\mu$m)")
        ax.margins(x=0)

        axcolor = 'lightgoldenrodyellow'
        axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)

        sfreq = Slider(axfreq, 'Light Offset', -self.L/2, +self.L/2, valinit=self.light_offset, valstep=self.dx,valfmt="%2.0e")


        def update(val):
            freq = np.float(sfreq.val)
            fd.create_source(offset=freq,plotOn=False)
            img.set_data(fd.calculate_propagation(plotOn=False))
            fig.canvas.draw_idle()


        sfreq.on_changed(update)

        plt.show()

        

    

if __name__ == "__main__":
    # %matplotlib qt
    fd=FdBpm()
    fd.NUM_SAMPLES=100
    fd.LENGTH=1.2E4
    fd.dy=1E-6
    fd.l_ambda=0.8E-6
    fd.L=60E-6
    fd.create_space()
    plotOn=True
    offset=13E-6
    fd.create_source(waist=5E-6,offset=-offset,plotOn=False)
    fd.create_guides(offset=-offset)
    fd.create_guides()
    fd.create_guides(offset=offset,plotOn=True)
    # fd.calculate_propagation()
    fd.plot_moving_source()