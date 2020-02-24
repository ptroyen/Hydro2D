import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

p = []
T = []
x = []
y = []
z = [] 


fparms = 'parameters.inp'
nc_col,nc_row = np.loadtxt(fparms)[:]


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


mx,my = np.meshgrid(cx,cy) 



import csv
with open('for_isolines0.9.csv') as csvfile:
  
    readCSV = csv.reader(csvfile, delimiter=',')
    headers = next(csvfile) 
    for row in readCSV:
        if float(row[7]) > 0.0:
            temp_p = row[0]
            temp_T = row[1]
            temp_x = row[5]
            temp_y = row[6]
            temp_z = row[7]
            p.append(float(temp_p))
            T.append(float(temp_T))
            x.append(float(temp_x))
            y.append(float(temp_y))
            z.append(float(temp_z))
        
        
        
# target grid to interpolate to
xi = np.arange(0.2,0.6,0.0001)
yi = np.arange(0.05,0.2,0.0001)
xi,yi = np.meshgrid(xi,yi)


# interpolate
#zi = griddata((x,y),z,(xi,yi),method='linear')
pi = griddata((x,y),p,(xi,yi),method='linear',fill_value=0.0)
Ti = griddata((x,y),T,(xi,yi),method='linear',fill_value=0.0)  

# np.savetxt("cubic.csv", (pi.flatten(), Ti.flatten()))

###PLOTS
# fig, ax = plt.subplots()
# CS = ax.contour(xi, yi, pi)
# ax.clabel(CS, inline=0, fontsize=10)
# ax.set_title('Simplest default with labels')

# tt = np.around(Ti, decimals=0, out=None)
# tt = tt.astype(int)

fig, ax = plt.subplots()
fig.set_size_inches(15, 15)
CS = ax.contour(xi, yi, Ti)
# CS = ax.contour(xi, yi, Ti,levels = [200,300,400,550,700,800,900,950,1000,1100,1150,1200],cmap= 'gray_r' )
# CS = ax.contourf(xi, yi, Ti,50)
fig.colorbar(CS, ax=ax, shrink=1)
#ax.set_ylim([0.07,0.17])
#ax.set_xlim([0.2,0.6])
#ax.clabel(CS, inline=1, fontsize=12)
ax.set_title('for temperature; cubic interpolation')

fig, ax = plt.subplots()
fig.set_size_inches(15, 15)
CS = ax.contour(xi, yi, Ti,6,cmap= 'gray_r' )
# CS = ax.contour(xi, yi, Ti,levels = [200,300,400,550,700,800,900,950,1000,1100,1150,1200],cmap= 'gray_r' )
# CS = ax.contourf(xi, yi, Ti,50)
#ax.set_ylim([0.07,0.17])
#ax.set_xlim([0.2,0.6])
fig.colorbar(CS, ax=ax, shrink=1)
# ax.clabel(CS, inline=1, fontsize=12)
ax.set_title('for temperature; cubic interpolation, 6 levels')

# tt = np.around(Ti, decimals=0, out=None)
# tt = tt.astype(int)
fig, ax = plt.subplots()
fig.set_size_inches(15, 15)
# CS = ax.contour(xi, yi, Ti,6,colors='k')
CS = ax.contour(xi, yi, Ti,levels = [200,300,550,700,800,900,950,1100,1150,1200],cmap='gray_r')
# CS = ax.contourf(xi, yi, Ti,50)
#ax.set_ylim([0.07,0.17])
#ax.set_xlim([0.2,0.6])
fig.colorbar(CS, ax=ax, shrink=1)
# ax.clabel(CS, inline=1, fmt= '%d',fontsize=9)
ax.set_title('for temperature; cubic interpolation, 10 levels specified')

fig,ax = plt.subplots()
fig.set_size_inches(15, 15)
plt.scatter(x, y, c=T , cmap='plasma')
ax.set_ylim([0.03,0.24])
ax.set_xlim([0.2,0.6])
fig.set_size_inches(15, 15)
CS = ax.contour(xi, yi, Ti,levels = [200,300,550,700,800,900,950,1100,1150,1200],colors='black')
#CS = ax.contourf(xi, yi, pi,150)
ax.set_ylim([0.03,0.24])
ax.set_xlim([0.2,0.6])
ax.set_title('for temperature; imported data and isolines compared')
plt.show()