import numpy as np
import matplotlib.pyplot as pyplot
import argparse


fname = "check.txt"




fnameout = 'out_check.txt'


with open(fname, 'r') as fin:
    data = fin.read().splitlines(True)
with open(fnameout, 'w') as fout:
    fout.writelines(data[1:])

header = data[0]
print(header)



X = np.loadtxt(fnameout)[:, 0]
Y = np.loadtxt(fnameout)[:, 1]
R = np.loadtxt(fnameout)[:, 2]
U = np.loadtxt(fnameout)[:, 3]
V = np.loadtxt(fnameout)[:, 4]
P = np.loadtxt(fnameout)[:, 5]
FLX = np.loadtxt(fnameout)[:, 6]
FRX = np.loadtxt(fnameout)[:, 7]
FLY = np.loadtxt(fnameout)[:, 8]
FRY = np.loadtxt(fnameout)[:, 9]


    
fparms = 'parameters.inp'
nc_col,nc_row = np.loadtxt(fparms)[:]

nc_col = int(nc_col) + 4
nc_row = int(nc_row) + 4

# ncx = nc_col
# ncy = nc_row

ncx = nc_col
ncy = nc_row

n_col = nc_col + 1
n_row = nc_row + 1

cx = np.zeros(nc_col)
cy = np.zeros(nc_row)

for j in range(nc_col):
    cx[j] = X[j]
    
for i in range(nc_row):
    cy[i] = Y[i*nc_col]


r = np.zeros([nc_row,nc_col])
u = np.zeros([nc_row,nc_col])
v = np.zeros([nc_row,nc_col])
p = np.zeros([nc_row,nc_col])
flx = np.zeros([nc_row,nc_col])
frx = np.zeros([nc_row,nc_col])
fly = np.zeros([nc_row,nc_col])
fry = np.zeros([nc_row,nc_col])

S = np.sqrt(U*U + V*V) ## caclulate speed as well







for i in range(nc_row):
    for j in range(nc_col):
        r[i][j] = R[i*nc_col+j]
        u[i][j] = U[i*nc_col+j]
        v[i][j] = V[i*nc_col+j]
        p[i][j] = P[i*nc_col+j]

        flx[i][j] = FLX[i*nc_col+j]
        frx[i][j] = FRX[i*nc_col+j]
        fly[i][j] = FLY[i*nc_col+j]
        fry[i][j] = FRY[i*nc_col+j]

        
mx,my = np.meshgrid(cx,cy)     


pyplot.figure(5)
pyplot.contourf(cx,cy,r,160,cmap = 'RdGy')
pyplot.colorbar()
# pyplot.legend()
pyplot.title("density")
pyplot.show()

pyplot.figure(6)
pyplot.contourf(cx,cy,p,60,cmap = 'RdGy')
pyplot.colorbar()
# pyplot.legend()
pyplot.title("pressure")
pyplot.show()

pyplot.figure(7)
pyplot.contourf(cx,cy,u,80,cmap = 'RdGy')
pyplot.colorbar()
# pyplot.legend()
pyplot.title("velocity_x")
pyplot.show()

pyplot.figure(8)
pyplot.contourf(cx,cy,v,80,cmap = 'RdGy')
pyplot.colorbar()
# pyplot.legend()
pyplot.title("velocity_y")
pyplot.show()

pyplot.figure(99)
pyplot.contour(cx,cy,flx,80,cmap = 'RdGy')
# pyplot.scatter(cx,cy,flx)
# pyplot.colorbar()
# pyplot.legend()
pyplot.title("flx")
pyplot.show()

pyplot.figure(100)
pyplot.contourf(cx,cy,frx,160,cmap = 'RdGy')
pyplot.colorbar()
# pyplot.legend()
pyplot.title("frx")
pyplot.show()

pyplot.figure(54)
pyplot.contourf(cx,cy,fly,80,cmap = 'RdGy')
pyplot.colorbar()
# pyplot.legend()
pyplot.title("fly")
pyplot.show()

pyplot.figure(56)
pyplot.contourf(cx,cy,fry,160,cmap = 'RdGy')
pyplot.colorbar()
# pyplot.legend()
pyplot.title("fry")
pyplot.show()
