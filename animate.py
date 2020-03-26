import numpy as np
import matplotlib.pyplot as pyplot

fparms = 'parameters.inp'
nc_col,nc_row = np.loadtxt(fparms)[:]

nc_col = int(nc_col)
nc_row = int(nc_row)

ncx = nc_col
ncy = nc_row

n_col = nc_col + 1
n_row = nc_row + 1

cx = np.zeros(nc_col)
cy = np.zeros(nc_row)

X = np.loadtxt("testout.txt")[:, 0]
Y = np.loadtxt("testout.txt")[:, 1]

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

# f = open(fname,'r')
# lines = f.readlines()[1:]
# f.close()

# with open(fname) as f:
#     next(f)
#     for line in f:
#         #do something

i = 50
while i < 300000 :

    nam = str(i)
    fname = nam+".txt"

    print(fname)


    with open(fname, 'r') as fin:
        data = fin.read().splitlines(True)
    with open(fname, 'w') as fout:
        fout.writelines(data[1:])

    i = i+ 50

i = 50
while i < 300000 :

    nam = str(i)
    fname = nam+".txt"

    print(fname)

    R = np.loadtxt(fname)[:, 2]
    U = np.loadtxt(fname)[:, 3]
    V = np.loadtxt(fname)[:, 4]
    P = np.loadtxt(fname)[:, 5]
    E = np.loadtxt(fname)[:, 6]



    for k in range(nc_row):
        for l in range(nc_col):
            r[i][j] = R[i*nc_col+j]
            u[i][j] = U[i*nc_col+j]
            v[i][j] = V[i*nc_col+j]
            p[i][j] = P[i*nc_col+j]
            e[i][j] = E[i*nc_col+j]
            

    mx,my = np.meshgrid(cx,cy)  

    pyplot.figure(5)
    pyplot.contourf(cx,cy,r,160,cmap = 'RdGy')
    pyplot.colorbar()
    # pyplot.legend()
    pyplot.title("density")
    pyplot.draw()

    pyplot.figure(7)
    pyplot.contourf(cx,cy,u,80,cmap = 'RdGy')
    pyplot.colorbar()
    # pyplot.legend()s
    pyplot.title("velocity_x")
    pyplot.draw()

    i = i + 50    


    
# import numpy as np
# from matplotlib import pyplot as plt
# from matplotlib.animation import FuncAnimation
# plt.style.use('seaborn-pastel')


# fig = plt.figure()
# ax = plt.axes(xlim=(0, 4), ylim=(-2, 2))
# line, = ax.plot([], [], lw=3)

# def init():
#     line.set_data([], [])
#     return line,
# def animate(i):
#     x = np.linspace(0, 4, 1000)
#     y = np.sin(2 * np.pi * (x - 0.01 * i))
#     line.set_data(x, y)
#     return line,

# anim = FuncAnimation(fig, animate, init_func=init,
#                                frames=200, interval=20, blit=True)


# anim.save('sine_wave.gif', writer='imagemagick')


        

   

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






# pyplot.figure(55)
# pyplot.contour(cx,cy,r,10)#,cmap = 'RdGy')
# # pyplot.colorbar()
# # pyplot.legend()
# pyplot.title("density")
# pyplot.show()

# pyplot.figure(6)
# pyplot.contourf(cx,cy,p,60,cmap = 'RdGy')
# pyplot.colorbar()
# # pyplot.legend()
# pyplot.title("pressure")
# pyplot.show()


# pyplot.figure(7)
# pyplot.contourf(cx,cy,u,80,cmap = 'RdGy')
# pyplot.colorbar()
# # pyplot.legend()
# pyplot.title("velocity_x")
# pyplot.show()


# pyplot.figure(8)
# pyplot.contourf(cx,cy,v,80,cmap = 'RdGy')
# pyplot.colorbar()
# # pyplot.legend()
# pyplot.title("velocity_y")
# pyplot.show()



# pyplot.figure(99)
# pyplot.contourf(cx,cy,s,80,cmap = 'RdGy')
# pyplot.colorbar()
# # pyplot.legend()
# pyplot.title("Speed")
# pyplot.show()

# pyplot.figure(100)
# pyplot.contourf(cx,cy,e,80,cmap = 'RdGy')
# pyplot.colorbar()
# # pyplot.legend()
# pyplot.title("spc_energy")
# pyplot.show()



# # pyplot.figure(9)
        
# # pyplot.plot(cx,r[int(ncx/2),:],label='density along x')
# # pyplot.legend()
# # pyplot.show()


# # pyplot.figure(10)
# # pyplot.plot(cy,v[:,int(ncy/2)],label='velocity_y cut along y')
# # pyplot.legend()
# # pyplot.show()

# # pyplot.figure(12)
        
# # pyplot.plot(cy,r[:,int(ncy/2)],label='density along y')
# # pyplot.legend()
# # pyplot.show()


# # pyplot.figure(13)
# # pyplot.plot(cx,v[int(ncx/2),:],label='velocity_y cut along x')
# # pyplot.legend()
# # pyplot.show()

# # pyplot.figure(14)
# # pyplot.plot(cx,u[int(ncx/2),:],label='velocity_x cut along x')
# # pyplot.legend()
# # pyplot.show()

# # pyplot.figure(142)
# # pyplot.plot(cx,p[int(ncx/2),:],label='pressure cut along x')
# # pyplot.legend()
# # pyplot.show()
        

        
# # pyplot.figure(15)
# # pyplot.contourf(cx,cy,speed,80,cmap = 'RdGy')
# # pyplot.colorbar()
# # # pyplot.legend()
# # pyplot.title("velocity_r")
# # pyplot.show()
        

# # flat_x = numpy.ndarray.flatten(mesh_x)
# # flat_y = numpy.ndarray.flatten(mesh_y)
# # flat_r = numpy.ndarray.flatten(r)
# # flat_u = numpy.ndarray.flatten(u)
# # flat_v = numpy.ndarray.flatten(v)
# # flat_p = numpy.ndarray.flatten(p)


# # dump = numpy.stack((flat_x,flat_y,flat_r,flat_u,flat_v,flat_p),axis=-1)
# # numpy.savetxt('output2D.txt',dump)  


# # pyplot.figure(16)
        
# # pyplot.plot(cy,v[:,int(ncy/2)],label='$V_y$ along x = 0.9')
# # pyplot.plot(cx,u[int(ncx/2),:],label='$V_x$ along y = 0.9')
# # pyplot.legend()
# # pyplot.show()

# # pyplot.figure(17)

# # pyplot.plot(cx,r[int(ncx/2),:],label='density along y = 0.9')
# # pyplot.plot(cy,r[:,int(ncy/2)],label='density along x = 0.9')
# # pyplot.legend()
# # pyplot.show()

# # pyplot.figure(18)

# # pyplot.plot(cx,e[int(ncx/2),:],label='energy along y = 0.9')
# # pyplot.plot(cy,e[:,int(ncy/2)],label='energy along x = 0.9')
# # pyplot.legend()
# # pyplot.show()


