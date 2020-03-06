import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import style

global i


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
k=20000
# fig, (ax1,ax2) = plt.subplots(1,2)



plt.ion()
plt.show()

# def animate(a):
while k < 600000 :
        fig = plt.figure(1)
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)

        ax1.title.set_text('Density')
        ax2.title.set_text('Velocity_X')
        
        nam = str(k)
        fname = nam+".txt"

        print(fname)


        with open(fname, 'r') as fin:
            data = fin.read().splitlines(True)
        with open('out_'+fname, 'w') as fout:
            fout.writelines(data[1:])


        # print(a)
        print("************")


        fname = 'out_'+fname

        # print(fname)

        R = np.loadtxt(fname)[:, 2]
        U = np.loadtxt(fname)[:, 3]
        # V = np.loadtxt(fname)[:, 4]
        # P = np.loadtxt(fname)[:, 5]
        # E = np.loadtxt(fname)[:, 6]



        for i in range(nc_row):
            for j in range(nc_col):
                r[i][j] = R[i*nc_col+j]
                u[i][j] = U[i*nc_col+j]
                # v[i][j] = V[i*nc_col+j]
                # p[i][j] = P[i*nc_col+j]
                # e[i][j] = E[i*nc_col+j]
                

        # mx,my = np.meshgrid(cx,cy)  

        # pyplot.figure(5)
        
        
        # pyplot.legend()
       
        ax1.clear()
        
        p1 = ax1.contourf(cx,cy,r,160,cmap = 'RdGy')
        cb1=pyplot.colorbar(p1,ax=ax1)
        # cb1.remove()
        # pyplot.draw()
        

        # pyplot.figure(7)
        
        
        # pyplot.legend()s
        ax2.clear()
        
        p2 = ax2.contourf(cx,cy,u,80,cmap = 'RdGy')
        cb2=pyplot.colorbar(p2,ax=ax2)
        # cb2.remove()
        pyplot.draw()
        plt.pause(0.0001)
        k = k+ 500

        fig.clear()
        

        



# f = open(fname,'r')
# lines = f.readlines()[1:]
# f.close()

# with open(fname) as f:
#     next(f)
#     for line in f:
#         #do somethin

# (ax1,ax2) = fig.add_subplot(1,2,2)

        

# def animate(i):
#     graph_data = open('example.txt','r').read()
#     lines = graph_data.split('\n')
#     xs = []
#     ys = []
#     for line in lines:
#         if len(line) > 1:
#             x, y = line.split(',')
#             xs.append(float(x))
#             ys.append(float(y))
#     ax1.clear()
#     ax1.plot(xs, ys)



# anim = animation.FuncAnimation(fig, animate,frames=100, interval=20, blit=True)
# plt.show()




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


