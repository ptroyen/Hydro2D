import numpy as np
import matplotlib.pyplot as pyplot
import argparse


dstart = 500

parser = argparse.ArgumentParser(description='Process Plot File')
parser.add_argument('name', metavar='N', type=int,default=dstart, help = 'Start Numerical Value of the Filename for Plotting')
parser.add_argument('locate', metavar='L', type=int,default=0, help = 'Plot Initial or Final [-1 for initial 1 for final]')
# parser.add_argument('incre', metavar='I', type=int,default=dincre, help = 'Stride in Numerical Value for Increment')
# parser.add_argument('add', metavar='A', type=int,default=0, help = 'Animate Initial Addition')

inp = parser.parse_args()

X = np.loadtxt("testout2.txt")[:, 0]
Y = np.loadtxt("testout2.txt")[:, 1]

loc = inp.locate
nam = inp.name
sname = str(nam)
fname = sname+'.txt'
fnameout = 'out_plot.txt'

if ( loc == 1):
    
    R = np.loadtxt("testout.txt")[:, 2]
    U = np.loadtxt("testout.txt")[:, 3]
    V = np.loadtxt("testout.txt")[:, 4]
    P = np.loadtxt("testout.txt")[:, 5]
    E = np.loadtxt("testout.txt")[:, 6]
elif(loc == -1):
    fname = 'testout2.txt'
    R = np.loadtxt(fname)[:, 2]
    U = np.loadtxt(fname)[:, 3]
    V = np.loadtxt(fname)[:, 4]
    P = np.loadtxt(fname)[:, 5]
    E = np.loadtxt(fname)[:, 6]
elif(loc == 2):
    with open(fname, 'r') as fin:
        data = fin.read().splitlines(True)
    with open(fnameout, 'w') as fout:
        fout.writelines(data[1:])

 
    

    R = np.loadtxt(fnameout)[:, 2]
    U = np.loadtxt(fnameout)[:, 3]
    V = np.loadtxt(fnameout)[:, 4]
    P = np.loadtxt(fnameout)[:, 5]
    E = np.loadtxt(fnameout)[:, 6]

else:
    with open('add'+fname, 'r') as fin:
        data = fin.read().splitlines(True)
    with open(fnameout, 'w') as fout:
        fout.writelines(data[1:])

    R = np.loadtxt(fnameout)[:, 2]
    U = np.loadtxt(fnameout)[:, 3]
    V = np.loadtxt(fnameout)[:, 4]
    P = np.loadtxt(fnameout)[:, 5]
    E = np.loadtxt(fnameout)[:, 6]

    # X = np.loadtxt("initial.txt")[:, 0]
    # Y = np.loadtxt("initial.txt")[:, 1]
    # R = np.loadtxt("initial.txt")[:, 2]
    # U = np.loadtxt("initial.txt")[:, 3]
    # V = np.loadtxt("initial.txt")[:, 4]
    # P = np.loadtxt("initial.txt")[:, 5]
    # E = np.loadtxt("initial.txt")[:, 6]

    # X = np.loadtxt("testout2.txt")[:, 0]
    # Y = np.loadtxt("testout2.txt")[:, 1]
    # R = np.loadtxt("testout2.txt")[:, 2]
    # U = np.loadtxt("testout2.txt")[:, 3]
    # V = np.loadtxt("testout2.txt")[:, 4]
    # P = np.loadtxt("testout2.txt")[:, 5]
    # E = np.loadtxt("testout2.txt")[:, 6]

# # Time =  0.0000068211804839
#     X = np.loadtxt("4200.txt")[:, 0]
#     Y = np.loadtxt("4200.txt")[:, 1]
#     R = np.loadtxt("4200.txt")[:, 2]
#     U = np.loadtxt("4200.txt")[:, 3]
#     V = np.loadtxt("4200.txt")[:, 4]
#     P = np.loadtxt("4200.txt")[:, 5]
#     E = np.loadtxt("4200.txt")[:, 6]




# X = np.loadtxt("2.txt")[:, 0]
# Y = np.loadtxt("2.txt")[:, 1]
# R = np.loadtxt("2.txt")[:, 2]
# U = np.loadtxt("2.txt")[:, 3]
# V = np.loadtxt("2.txt")[:, 4]
# P = np.loadtxt("2.txt")[:, 5]
# E = np.loadtxt("2.txt")[:, 6]

## COL 4
## ROW 3

Time = data[0]
print(Time)
    
fparms = 'parameters.inp'
nc_col,nc_row = np.loadtxt(fparms)[:]

nc_col = int(nc_col)
nc_row = int(nc_row)

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
e = np.zeros([nc_row,nc_col])
s = np.zeros([nc_row,nc_col])
S = np.sqrt(U*U + V*V) ## caclulate speed as well







for i in range(nc_row):
    for j in range(nc_col):
        r[i][j] = R[i*nc_col+j]
        u[i][j] = U[i*nc_col+j]
        v[i][j] = V[i*nc_col+j]
        p[i][j] = P[i*nc_col+j]
        e[i][j] = E[i*nc_col+j]
        s[i][j] = S[i*nc_col+j]
        
    
    
    
    ## PLOTS HERE
        

p_rho = p/r ## some sort of energy
# print(np.min(r))

mx,my = np.meshgrid(cx,cy)     

## plot mesh
# pyplot.figure(211)

# pyplot.plot(cx,mx,'k')
# pyplot.plot(my.T,cy,'k')
# ## pyplot.plot(cy)
# pyplot.legend()
# pyplot.show()

# pyplot.figure()
# pyplot.plot(mx,my,'b')
# pyplot.plot(mx.T,my.T,'b')
# pyplot.title("Algebraic Grid")
# pyplot.show()


levels = np.linspace(250000,260000,200)
pyplot.figure(100)
pyplot.contourf(cx,cy,e,levels = levels,cmap = 'plasma')
pyplot.colorbar()
# pyplot.legend()
pyplot.title("spc_energy")
pyplot.show()

pyplot.figure(5)
pyplot.contourf(cx,cy,r,160,cmap = 'RdGy')
pyplot.colorbar()
# pyplot.legend()
pyplot.title("density")
pyplot.show()



pyplot.figure(55)
pyplot.contour(cx,cy,r,10)#,cmap = 'RdGy')
# pyplot.colorbar()
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
pyplot.contourf(cx,cy,s,80,cmap = 'RdGy')
pyplot.colorbar()
# pyplot.legend()
pyplot.title("Speed")
pyplot.show()






# pyplot.figure(9)
        
# pyplot.plot(cx,r[int(ncx/2),:],label='density along x')
# pyplot.legend()
# pyplot.show()


# pyplot.figure(10)
# pyplot.plot(cy,v[:,int(ncy/2)],label='velocity_y cut along y')
# pyplot.legend()
# pyplot.show()

# pyplot.figure(12)
        
# pyplot.plot(cy,r[:,int(ncy/2)],label='density along y')
# pyplot.legend()
# pyplot.show()


# pyplot.figure(13)
# pyplot.plot(cx,v[int(ncx/2),:],label='velocity_y cut along x')
# pyplot.legend()
# pyplot.show()

# pyplot.figure(14)
# pyplot.plot(cx,u[int(ncx/2),:],label='velocity_x cut along x')
# pyplot.legend()
# pyplot.show()

# pyplot.figure(142)
# pyplot.plot(cx,p[int(ncx/2),:],label='pressure cut along x')
# pyplot.legend()
# pyplot.show()
        

        
# pyplot.figure(15)
# pyplot.contourf(cx,cy,speed,80,cmap = 'RdGy')
# pyplot.colorbar()
# # pyplot.legend()
# pyplot.title("velocity_r")
# pyplot.show()
        

# flat_x = numpy.ndarray.flatten(mesh_x)
# flat_y = numpy.ndarray.flatten(mesh_y)
# flat_r = numpy.ndarray.flatten(r)
# flat_u = numpy.ndarray.flatten(u)
# flat_v = numpy.ndarray.flatten(v)
# flat_p = numpy.ndarray.flatten(p)


# dump = numpy.stack((flat_x,flat_y,flat_r,flat_u,flat_v,flat_p),axis=-1)
# numpy.savetxt('output2D.txt',dump)  


# pyplot.figure(16)
        
# pyplot.plot(cy,v[:,int(ncy/2)],label='$V_y$ along x = 0.9')
# pyplot.plot(cx,u[int(ncx/2),:],label='$V_x$ along y = 0.9')
# pyplot.legend()
# pyplot.show()

# pyplot.figure(17)

# pyplot.plot(cx,r[int(ncx/2),:],label='density along y = 0.9')
# pyplot.plot(cy,r[:,int(ncy/2)],label='density along x = 0.9')
# pyplot.legend()
# pyplot.show()

# pyplot.figure(18)

# pyplot.plot(cx,e[int(ncx/2),:],label='energy along y = 0.9')
# pyplot.plot(cy,e[:,int(ncy/2)],label='energy along x = 0.9')
# pyplot.legend()
# pyplot.show()


