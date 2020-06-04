#%%
import numpy as np
import time as t
import matplotlib.pyplot as plt

x1 = -60e-6
x2 = 60e-6
num_samples = 50
dx = (x2-x1)/num_samples
dy = 1e-6
x = np.linspace(x1, x2-dx, num_samples)
W0 = 8e-6
n_index = np.ones((num_samples,))
nmax = 1
nmin = 1
n_bar = (nmax+nmin)/2
l_ambda = 0.8e-6
k0 = 2*np.pi/l_ambda

#%% Gauss
modo = np.array(np.exp(-(x/W0)**2))
modo2=modo

k = k0*n_index
k_bar = k*n_bar

length=1000

h = dy
ro = dy/(dx**2)
A = 1j/(2*k_bar)
B = 1j*(k**2-k_bar**2)/(2*k_bar)
a = -ro*A
b = 2*(1+ro*A)-h*B
c = a

matrix = np.zeros((num_samples,num_samples),dtype='complex_')
# matrix_original = np.zeros((num_samples,num_samples),dtype='complex_')

start=t.time()

#%% Tridiagonal Matrix

index=np.arange(0,num_samples-1)
matrix[index,index+1]=a[index]

index=np.arange(1,num_samples)
matrix[index,index-1]=c[index]

index=np.arange(num_samples)
matrix[index,index]=b[index]

# print((matrix_original==matrix).all())
#%% Propagation
modo = np.array(np.exp(-(x/W0)**2))
zz = np.zeros((length,num_samples))
d = np.zeros((num_samples),dtype='complex_')


for n in range(length):         # cannot be transferred to numpy cause of np.linalg.solve in between lines
    index=np.arange(1,num_samples-1)
    d[index]=(2*(1-ro*A[index])+h*B[index])*modo[index]+ro*A[index]*(modo[index-1]+modo[index+1])
    d[0] = (2*(1-ro*A[0])+h*B[0])*modo[0]+ro*A[0]*(modo[1])
    d[-1] = (2*(1-ro*A[-1])+h*B[-1])*modo[-1]+ro*A[-1]*(modo[-2])

    # modo = np.linalg.lstsq(matrix,d)[0]
    # modo = solveit(matrix,d)
    modo = np.linalg.solve(matrix,d) #fastest option
    zz[n,:] = np.abs(modo)

#%%

# print(modo[0])

end = t.time()
print("Time elapsed:%0.4f seconds"%(end-start))

#%%

zz=zz[::np.int(np.floor(length/num_samples)),:]

plt.imshow(zz,cmap="jet",interpolation='bilinear')
plt.show()

# %%
