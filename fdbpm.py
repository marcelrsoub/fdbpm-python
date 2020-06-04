#%%
import numpy as np
import time as t
import matplotlib.pyplot as plt

x1 = -60e-6
x2 = 60e-6
num_samples = 101
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

#%%Gauss
modo = np.array(np.exp(-(x/W0)**2))

k = k0*n_index
k_bar = k*n_bar

length=1000

#matrices
# xx,yy = np.meshgrid(np.arange(x1,x2,dx),np.arange(dy,1000e-6,dy))
zz = np.zeros((length,num_samples))
h = dy
ro = dy/(dx**2)
A = 1j/(2*k_bar)
B = 1j*(k**2-k_bar**2)/(2*k_bar)
a = -ro*A
b = 2*(1+ro*A)-h*B
c = a
d = np.zeros((num_samples),dtype='complex_')
matrix = np.zeros((num_samples,num_samples),dtype='complex_')
# matrix_original = np.zeros((num_samples,num_samples),dtype='complex_')

start=t.time()

#%%
# % --------- Generacion de la matriz tridiagonal ---------

for m in range(num_samples):
    if ((m>0) and (m<num_samples-1)):
        matrix_original[m,m-1] = a[m]
        matrix_original[m,m] = b[m]
        matrix_original[m,m+1] = c[m]
    else:
        matrix_original[0,0] = b[0]
        matrix_original[0,1] = c[0]
        matrix_original[num_samples-1,num_samples-2] = a[num_samples-1]
        matrix_original[num_samples-1,num_samples-1] = b[num_samples-1]
matrix=matrix_original
# index=np.arange(0,num_samples-1)
# matrix[index,index+1]=a[index]

# index=np.arange(1,num_samples)
# matrix[index,index-1]=c[index]

# index=np.arange(num_samples)
# matrix[index,index]=b[index]

# print((matrix_original==matrix).all())
#%%# % --------- Ciclo Principal de Propagacion ---------  
from itertools import combinations

def solveit(A,b):
    num_vars = A.shape[1]
    rank = np.linalg.matrix_rank(A)
    if rank == num_vars:              
        sol = np.linalg.lstsq(A, b)[0]    # not under-determined
    else:
        for nz in combinations(range(num_vars), rank):    # the variables not set to zero
            try: 
                sol = np.zeros((num_vars, 1))  
                sol[nz, :] = np.asarray(np.linalg.solve(A[:, nz], b))
                print(sol)
            except np.linalg.LinAlgError:     
                pass                    # picked bad variables, can't solve
    return sol
for n in range(length):
    for m in range(num_samples):
        if ((m>0) and (m<num_samples-1)):
            d[m] = (2*(1-ro*A[m])+h*B[m])*modo[m]+ro*A[m]*(modo[m-1]+modo[m+1])
        else:
            d[0] = (2*(1-ro*A[0])+h*B[0])*modo[0]+ro*A[0]*(modo[1])
            d[-1] = (2*(1-ro*A[-1])+h*B[-1])*modo[-1]+ro*A[-1]*(modo[-2])

    # modo = np.linalg.lstsq(matrix,d)[0]
    # modo = solveit(matrix,d)
    modo = np.linalg.solve(matrix,d)
    # zz[n,:] = modo
    zz[n,:] = np.abs(modo)
    # zz=np.abs(zz)

# print(modo[0])

end = t.time()
print("Time elapsed:%0.4f seconds"%(end-start))

#%%
num_vars = matrix.shape[1]
rank = np.linalg.matrix_rank(matrix)
if rank == num_vars:
    print("equal")

zz=zz[::10,:]
if (zz[0]==zz[-1]).all():
    print("Shit, wrong")
else:
    print("Hhhmmmm")
plt.imshow(zz,cmap="jet",interpolation='bilinear')
plt.show()

# %%
