import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D # <--- This is important for 3d plotting 
from scipy.interpolate import griddata

PI = np.pi

### MESH PART
# Mesh input parameters
fparms = 'parameters.inp'
nc_col,nc_row = np.loadtxt(fparms)[:]

nc_col = int(nc_col)
nc_row = int(nc_row)

ncx = nc_col
ncy = nc_row

cx = np.zeros(nc_col)
cy = np.zeros(nc_row)

X = np.loadtxt("mesh.txt")[:, 0]
Y = np.loadtxt("mesh.txt")[:, 1]

## CENTROIDS
for j in range(nc_col):
    cx[j] = X[j]
    
for i in range(nc_row):
    cy[i] = Y[i*nc_col]


# interpolate on this grid
xi,yi = np.meshgrid(cx,cy)

# The Dual Pulse Laser :: SECOND PULSE
w0 = 200.0e-6 
lam = 1064.0e-9 
f = 300.0e-3
Rl_rn = 0.1*PI*(w0*w0) / lam    ## Scale Down the length in this direction

xmid = cx[int(nc_col/2)]
ymid = cy[int(nc_row/2)]

#normalized one
x_nor_i = (xi-xmid) / Rl_rn
y_nor_i = (yi - ymid) / w0

## INTERPOLATION VALUES EXTRACT
# Import The Files
## ---------------------------
pre1 = 'I_Pulse2_'
pre2 = 'Rho_e_'
ext = '.txt'
p3d = '.p3d'

## LOOP for diffrent files
for i in range(8):

    index = i*5
    str_i = str(index)

    fname1 = pre1+str_i
    fname2 = pre2+str_i

    fout = 'out_'+str_i + ext
    # fout2 = 'out_'+fname2 + ext


    Z_r = np.loadtxt(fname1+ext)[0,1:]
    R_r = np.loadtxt(fname1+ext)[1:,0]
    In = np.loadtxt(fname1+ext)[1:,1:] #Value read, already in matrix form
    Ne = np.loadtxt(fname2+ext)[1:,1:]

    In_stack = np.vstack((np.flipud(In),In))
    Ne_stack = np.vstack((np.flipud(Ne),Ne))
    
    In_stack = 1.0e2 * In_stack                 ## Intensity value scaled up
    Ne_stack = 1.0e2 * Ne_stack                 ## Ne value scaled up
    
    R_r_stack = np.hstack((np.flipud(-R_r),R_r))

    X_r, Y_r = np.meshgrid(Z_r, R_r_stack)
    V_r1 = In_stack 
    V_r2 = Ne_stack

    points = np.vstack((X_r.ravel(), Y_r.ravel()))
    points = points.T

    val1 = V_r1.ravel()
    val2 = V_r2.ravel()


    points_i = np.vstack((x_nor_i.ravel(), y_nor_i.ravel()))
    points_i = points_i.T
    

    # # interpolate
    V_i1 = griddata(points,val1,points_i,method='linear',fill_value=np.min(In_stack))
    V_i2 = griddata(points,val2,points_i,method='linear',fill_value=np.min(Ne_stack))

    # Save the Values in a file to be used by Hydro2D
    
    ## First Intensity :: Then Number Density
    result = np.vstack((V_i1, V_i2))
    np.savetxt(fout, result.T)


    # FOR PLOTS __________________________________________________________________________________________
    # twd_V1 = np.reshape(V_i1,(nc_row,nc_col))
    # twd_V2 = np.reshape(V_i2,(nc_row,nc_col))


    # fig = plt.figure()
    # ax = fig.gca(projection='3d')
    # ax.plot_surface(xi,yi,twd_V1, cmap='plasma')
    # # cset = ax.contourf(X, Y, Z-fit, zdir='z', offset=-4.0e12, cmap='plasma')
    # ax.set_title("From Function -- Not-Normalized")
    # #    ax.set_ylim(0.03,0.04)
    # plt.show()

    # fig = plt.figure()
    # ax = fig.gca(projection='3d')
    # ax.plot_surface(xi,yi,twd_V2, cmap='plasma')
    # # cset = ax.contourf(X, Y, Z-fit, zdir='z', offset=-4.0e12, cmap='plasma')
    # ax.set_title("From Function -- Not-Normalized")
    # # ax.set_zlim(0,np.max(Z)+2)
    # plt.show()

    # # Plot the 3D figure of the fitted function and the residuals.
    # fig = plt.figure()
    # ax = fig.gca(projection='3d')
    # ax.plot_surface(x_nor_i,y_nor_i,twd_V1, cmap='plasma')
    # # cset = ax.contourf(X, Y, Z-fit, zdir='z', offset=-4.0e12, cmap='plasma')
    # ax.set_title("From Function -- Normalized")
    # # ax.set_zlim(0,np.max(Z)+2)
    # plt.show()

    # # Plot the 3D figure of the fitted function and the residuals.
    # fig = plt.figure()
    # ax = fig.gca(projection='3d')
    # ax.plot_surface(x_nor_i,y_nor_i,twd_V2, cmap='plasma')
    # # cset = ax.contourf(X, Y, Z-fit, zdir='z', offset=-4.0e12, cmap='plasma')
    # ax.set_title("From Function -- Normalized")
    # # ax.set_zlim(0,np.max(Z)+2)
    # plt.show()


    # # Plot the 3D figure of the fitted function and the residuals.
    # fig = plt.figure()
    # ax = fig.gca(projection='3d')
    # ax.plot_surface(X_r,Y_r,V_r1, cmap='plasma')
    # # cset = ax.contourf(X, Y, Z-fit, zdir='z', offset=-4.0e12, cmap='plasma')
    # ax.set_title("Data: In")
    # # ax.set_zlim(0,np.max(Z)+2)
    # plt.show()

    # # Plot the 3D figure of the fitted function and the residuals.
    # fig = plt.figure()
    # ax = fig.gca(projection='3d')
    # ax.plot_surface(X_r,Y_r,V_r2, cmap='plasma')
    # # cset = ax.contourf(X, Y, Z-fit, zdir='z', offset=-4.0e12, cmap='plasma')
    # ax.set_title("Data : Ne")
    # # ax.set_zlim(0,np.max(Z)+2)
    # plt.show()



# ## Check by reading the data file that you just wrote
# # mesh is in xi , yi :: [][]
# In_int = np.loadtxt("out_0.txt")[:, 0]
# Ne_int = np.loadtxt("out_0.txt")[:, 1]

# twd_V1 = np.reshape(In_int,(nc_row,nc_col))
# twd_V2 = np.reshape(Ne_int,(nc_row,nc_col))

# fig = plt.figure()
# ax = fig.gca(projection='3d')
# ax.plot_surface(xi,yi,twd_V1, cmap='plasma')
# # cset = ax.contourf(X, Y, Z-fit, zdir='z', offset=-4.0e12, cmap='plasma')
# ax.set_title("From Interpolation In -- Not-Normalized")
# # ax.set_ylim(0.03,0.04)
# plt.show()

# fig = plt.figure()
# ax = fig.gca(projection='3d')
# ax.plot_surface(xi,yi,twd_V2, cmap='plasma')
# # cset = ax.contourf(X, Y, Z-fit, zdir='z', offset=-4.0e12, cmap='plasma')
# ax.set_title("From Interpolation Ne -- Not-Normalized")
# # ax.set_zlim(0,np.max(Z)+2)
# plt.show()

# fig = plt.figure()
# ax = fig.gca(projection='3d')
# ax.plot_surface(X_r,Y_r,V_r1, cmap='plasma')
# # cset = ax.contourf(X, Y, Z-fit, zdir='z', offset=-4.0e12, cmap='plasma')
# ax.set_title("Data")
# # ax.set_zlim(0,np.max(Z)+2)
# plt.show()


# fig = plt.figure()
# ax = fig.gca(projection='3d')
# ax.plot_surface(X_r,Y_r,V_r2, cmap='plasma')
# # cset = ax.contourf(X, Y, Z-fit, zdir='z', offset=-4.0e12, cmap='plasma')
# ax.set_title("Data")
# # ax.set_zlim(0,np.max(Z)+2)
# plt.show()



# plt.plot(points_i[:,0],points_i[:,1],V_i1)


# twd_V1 = np.reshape(V_i1,(nc_row,nc_col))
# twd_V2 = np.reshape(V_i2,(nc_row,nc_col))

# # Plot the 3D figure of the fitted function and the residuals.
# fig = plt.figure()
# ax = fig.gca(projection='3d')
# ax.plot_surface(x_nor_i,y_nor_i,twd_V1, cmap='plasma')
# # cset = ax.contourf(X, Y, Z-fit, zdir='z', offset=-4.0e12, cmap='plasma')
# ax.set_title("From Function -- Normalized")
# # ax.set_zlim(0,np.max(Z)+2)
# plt.show()

# fig = plt.figure()
# ax = fig.gca(projection='3d')
# ax.plot_surface(xi,yi,twd_V1, cmap='plasma')
# # cset = ax.contourf(X, Y, Z-fit, zdir='z', offset=-4.0e12, cmap='plasma')
# ax.set_title("From Function -- Not-Normalized")
# # ax.set_zlim(0,np.max(Z)+2)
# plt.show()


# fig = plt.figure()
# ax = fig.gca(projection='3d')
# ax.plot_surface(X_r,Y_r,V_r1, cmap='plasma')
# # cset = ax.contourf(X, Y, Z-fit, zdir='z', offset=-4.0e12, cmap='plasma')
# ax.set_title("Data")
# # ax.set_zlim(0,np.max(Z)+2)
# plt.show()


# # Plot the 3D figure of the fitted function and the residuals.
# fig = plt.figure()
# ax = fig.gca(projection='3d')
# ax.plot_surface(x_nor_i,y_nor_i,twd_V2, cmap='plasma')
# # cset = ax.contourf(X, Y, Z-fit, zdir='z', offset=-4.0e12, cmap='plasma')
# ax.set_title("From Function")
# # ax.set_zlim(0,np.max(Z)+2)
# plt.show()

# fig = plt.figure()
# ax = fig.gca(projection='3d')
# ax.plot_surface(X_r,Y_r,V_r2, cmap='plasma')
# # cset = ax.contourf(X, Y, Z-fit, zdir='z', offset=-4.0e12, cmap='plasma')
# ax.set_title("Data")
# # ax.set_zlim(0,np.max(Z)+2)
# plt.show()

# ##PLOTS ///////////////////////////////////// END ------------------------------------------------
