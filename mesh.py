import numpy as np
import matplotlib.pyplot as pyplot

i = 'mesh' 
fname = str(i)+'.txt'


X = np.loadtxt(fname)[:, 0]
Y = np.loadtxt(fname)[:, 1]
   

## COL 4
## ROW 3

nc_row = 200
nc_col = 200 

ncx = nc_col
ncy = nc_row


n_col = nc_col + 1
n_row = nc_row + 1

x = np.zeros(n_col)
y = np.zeros(n_row)

for j in range(n_col):
    x[j] = X[j]
    
for i in range(n_row):
    y[i] = Y[i*n_col]

  

mx,my = np.meshgrid(x,y)     

## plot mesh
pyplot.figure(211)

pyplot.plot(x,mx,'k')
pyplot.plot(my.T,y,'k')
## pyplot.plot(cy)
#pyplot.legend()
pyplot.show()
