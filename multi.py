import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import datetime

# Import The Files
## ---------------------------
pre1 = 'I_Pulse2_'
pre2 = 'Rho_e_'
ext = '.txt'


fname = pre1+"35"
fname2 = pre2+"35"


fout = 'out_'+fname + ext
fout2 = 'out_'+fname2 + ext


fmainout =   'out_fit_' + fname + '_' + fname2 + ext


fout_for_c = fmainout


Z_r = np.loadtxt(fname+ext)[0,1:]
R_r = np.loadtxt(fname+ext)[1:,0]
In = np.loadtxt(fname+ext)[1:,1:]

# In = In * 100

Ne = np.loadtxt(fname2+ext)[1:,1:]

Z_rm, R_rm = np.meshgrid(Z_r, R_r)




In_stack = np.vstack((np.flipud(In),In))
Ne_stack = np.vstack((np.flipud(Ne),Ne))
R_r_stack = np.vstack((np.flipud(-R_r),R_r))

Z_rwm, R_rwm = np.meshgrid(Z_r, R_r_stack)

# Product of Electron No density and Intensity
Pr = In_stack*Ne_stack 


## For convenience
X, Y = np.meshgrid(Z_r, R_r_stack)

## X and Y contain the full Meshgrid now
xmin = np.min(X)
xmax = np.max(X)
ymin = np.min(Y)
ymax = np.max(Y)


## -------------------------- WHAT TO FIT PROVIDE ARRAY IN Z

Z = Pr
# Z = Ne_stack
## ------------------------------------------------------
# LASER DATA EXTRACTED Z = Z(X,Y) ---------------------------------------------
# Plot the Read Data Once
fig = plt.figure(1)
ax = fig.gca(projection='3d')
ax.plot_surface(X, Y, In_stack, cmap='plasma')
ax.set_title(fname)
ax.set_zlim(0,np.max(In_stack)+2)
plt.show()

fig = plt.figure(2)
ax = fig.gca(projection='3d')
ax.plot_surface(X, Y, Ne_stack, cmap='plasma')
ax.set_title(fname2)
# ax.set_zlim(0,np.max(Ne_stack)+2)
plt.show()


fig = plt.figure(3)
ax = fig.gca(projection='3d')
ax.plot_surface(X, Y, Pr, cmap='plasma')
ax.set_title(fname+'*'+fname2)
# ax.set_zlim(0,np.max(Ne_stack)+2)
plt.show()
##---------------------------------------------------------------
guess_prm_e = [] 
guess_prm_l = [] 
guess_prm_g = []

# Initial guesses to the fit for all Gaussians : x0, y0, xalpha, yalpha, A
guess_prm_g = [ (-1.5, -0.0, 2.4, 0.4, 2.0e10),
               (0.02, 0.0, 2.6, 2.6, 6.0e20),
              (0.24, 0.070, 1.40, 0.80, 2.0e28),
              (-1.2, -0.0, 2.4, 0.4, 2.0e30)
               
             ]


# guess_prm_g = [(0.0, 0.0, 2.6, 2.6, 6.0e26),
#               (0.24, 0.070, 1.40, 0.80, 6.0e26),
#              (-1.0, 0.01, 2.4, 0.40, 6.0e26),
#                (-0.90, -0.0, 2.4, 0.4, 6.0e30),
#              ]

## For lorentzian :  a , b , c , d , f , g
# guess_prm_l = [(-2.0e+2,  9.1e+1, -2.0e+1, -1.29e+1 , -4.39169913e+10, -2.33478842e+10)] 

## For the Exponent Term : x0,y0,A,B
guess_prm_e = [] 
# guess_prm_e = [(10.0,2.0,1.0e-2,1.0e1)] 

# Global Definition
gaus_no = np.size(guess_prm_g)//5
loren_no = np.size(guess_prm_l)//6
expo_no = np.size(guess_prm_e)//4

# Flatten the initial guess parameter list.
# p0 = np.ravel(guess_prms)
p0 = [p for prms in guess_prm_g for p in prms]
l0 = [p for prms in guess_prm_l for p in prms]
e0 = [p for prms in guess_prm_e for p in prms]

allpars = np.hstack((p0 , l0 , e0))
# allpars = p0


# Our functions to fit
def gaussian(x, y, x0, y0, xalpha, yalpha, A):
    return A * np.exp( -((x-x0)**2/(2*xalpha**2)) -((y-y0)**2/(2*yalpha**2)))

def lorentzian(x_in, y_in , a , b , c , d , f , g): 
    temp = a / (1.0 + np.power((x_in-b)/c, 2.0)) + d * (1.0 + np.power((y_in-f)/g, 2.0))
    return temp

def exponent(x,y,x0,y0,A,B):
    ## Beer's Law for Intensity degeneraiton
    res = np.exp(-(A**2)*(y-y0)**2/(np.power((np.sqrt(1+((B*(x-x0))**2))) , 2)))
    return res
    
    



# This is the callable that is passed to curve_fit. M is a (2,N) array
# where N is the total number of data points in Z, which will be ravelled
# to one dimension.
def _fun_fit(M,*args):
    ## If you change the number of functions then change here too
    
    x, y = M
    arr = np.zeros(x.shape)
    for i in range(gaus_no): # // gives quotient
       arr += gaussian(x, y, *args[i*5:i*5+5]) # pointer used

    for m in range(loren_no): # // gives quotient :: Lorentzian requires 6 parameters
        # i = gaus_no + m 
        # print("Inisde lorentzian")
        # print(m)
        arr += lorentzian(x, y, *args[(gaus_no*5+m*6):((gaus_no*5+m*6)+6)]) # pointer used
        
    for n in range(expo_no): # // gives quotient
       arr = arr*exponent(x, y, *args[(loren_no*6+gaus_no*5+n*4):((loren_no*6+gaus_no*5+n*4)+4)]) # pointer used
    return arr





## 5 gaussian
'''
[-5.09458761e-01  4.36818178e-01  7.09277273e-01  3.59632979e-01
  1.99006865e+13  1.24063296e+00  4.58714805e-06  8.80736004e-01
  1.19576287e+00  6.58656642e+12 -2.05679718e+00 -2.97481466e-05
  6.59309236e-01  1.24932616e+00  6.10484320e+12 -5.09454333e-01
 -4.36823874e-01  7.09275158e-01  3.59629106e-01  1.99005693e+13]


RECORDED DATA FOR USE WITH 4 GAUSSIANS

 //// For I_Pulse2_0*Rho_e_0
gaus1.x0 = -2.470534e-01 ; gaus1.y0 = 2.465843e-07 ;  gaus1.sigma_x = 5.314265e-01 ; gaus1.sigma_y = 3.952746e-01 ;  gaus1.A = 3.215602e+25 ; 
gaus2.x0 = 7.488499e-01 ; gaus2.y0 = 1.407520e-06 ;  gaus2.sigma_x = 4.782765e-01 ; gaus2.sigma_y = 7.270592e-01 ;  gaus2.A = 8.904928e+24 ; 
gaus3.x0 = -1.557827e+00 ; gaus3.y0 = 1.342642e-05 ;  gaus3.sigma_x = 6.963636e-01 ; gaus3.sigma_y = 9.185818e-01 ;  gaus3.A = 7.962736e+24 ; 
gaus4.x0 = 1.913055e+00 ; gaus4.y0 = -1.765879e-05 ;  gaus4.sigma_x = -6.388076e-01 ; gaus4.sigma_y = -1.268513e+00 ;  gaus4.A = 4.262157e+24 ;


 //// For I_Pulse2_5*Rho_e_5
gaus1.x0 = -5.502246e+04 ; gaus1.y0 = -2.675539e+04 ;  gaus1.sigma_x = -1.023134e+06 ; gaus1.sigma_y = -2.805202e+03 ;  gaus1.A = -2.014646e+27 ; 
gaus2.x0 = -9.483050e+03 ; gaus2.y0 = -4.721620e+02 ;  gaus2.sigma_x = 1.498227e+03 ; gaus2.sigma_y = 5.811655e+02 ;  gaus2.A = 1.964729e+27 ; 
gaus3.x0 = 8.176103e-02 ; gaus3.y0 = -2.229833e+04 ;  gaus3.sigma_x = -2.908805e+00 ; gaus3.sigma_y = -7.419357e+03 ;  gaus3.A = 1.690062e+27 ; 
gaus4.x0 = -2.843474e-01 ; gaus4.y0 = -5.576152e-07 ;  gaus4.sigma_x = 4.052929e-01 ; gaus4.sigma_y = 1.234419e-01 ;  gaus4.A = 5.890770e+27 ; 

 //// For I_Pulse2_10*Rho_e_10
gaus1.x0 = -1.639137e+02 ; gaus1.y0 = -1.313185e+02 ;  gaus1.sigma_x = -6.487636e+02 ; gaus1.sigma_y = 8.851845e+02 ;  gaus1.A = 1.368640e+26 ; 
gaus2.x0 = -3.404239e-01 ; gaus2.y0 = -6.296805e-06 ;  gaus2.sigma_x = 4.144946e-01 ; gaus2.sigma_y = 1.173919e-01 ;  gaus2.A = 5.423196e+31 ; 
gaus3.x0 = -3.408081e-01 ; gaus3.y0 = -6.299272e-06 ;  gaus3.sigma_x = 4.089790e-01 ; gaus3.sigma_y = 1.152694e-01 ;  gaus3.A = -9.754017e+31 ; 
gaus4.x0 = -3.414428e-01 ; gaus4.y0 = -6.263239e-06 ;  gaus4.sigma_x = 4.019238e-01 ; gaus4.sigma_y = 1.127181e-01 ;  gaus4.A = 4.364736e+31 ; 

 //// For I_Pulse2_15*Rho_e_15   TIME::
gaus1.x0 = 1.687573e+05 ; gaus1.y0 = -5.344201e-02 ;  gaus1.sigma_x = 3.651564e+04 ; gaus1.sigma_y = -2.173010e+00 ;  gaus1.A = 2.555828e+31 ; 
gaus2.x0 = -3.819313e-01 ; gaus2.y0 = -1.575934e-09 ;  gaus2.sigma_x = 4.314179e-01 ; gaus2.sigma_y = -1.712695e-01 ;  gaus2.A = 1.155222e+30 ; 
gaus3.x0 = -3.566027e-01 ; gaus3.y0 = 1.357766e-06 ;  gaus3.sigma_x = 2.396164e-01 ; gaus3.sigma_y = 8.289540e-02 ;  gaus3.A = 7.931191e+29 ; 
gaus4.x0 = -5.109825e+02 ; gaus4.y0 = -5.336237e-01 ;  gaus4.sigma_x = -6.206727e+01 ; gaus4.sigma_y = -6.653129e+00 ;  gaus4.A = -2.463689e+31 ; 


 //// For I_Pulse2_20*Rho_e_20   TIME::
gaus1.x0 = -6.658634e+03 ; gaus1.y0 = -4.005392e+03 ;  gaus1.sigma_x = 1.473353e+05 ; gaus1.sigma_y = -1.813691e+03 ;  gaus1.A = 2.051803e+26 ; 
gaus2.x0 = 8.261236e+03 ; gaus2.y0 = 4.477484e-02 ;  gaus2.sigma_x = 7.413221e+03 ; gaus2.sigma_y = -4.464951e-01 ;  gaus2.A = 1.248243e+29 ; 
gaus3.x0 = -5.194304e+00 ; gaus3.y0 = 6.363464e-01 ;  gaus3.sigma_x = 3.652037e+01 ; gaus3.sigma_y = -1.285163e-01 ;  gaus3.A = -1.561502e+28 ; 
gaus4.x0 = -3.837622e-01 ; gaus4.y0 = -1.906526e-04 ;  gaus4.sigma_x = -3.217517e-01 ; gaus4.sigma_y = 1.457328e-01 ;  gaus4.A = 5.384729e+30 ; 

 //// For I_Pulse2_25*Rho_e_25   TIME::
gaus1.x0 = 1.000000e-01 ; gaus1.y0 = -0.000000e+00 ;  gaus1.sigma_x = 2.400000e+00 ; gaus1.sigma_y = 4.000000e-01 ;  gaus1.A = 2.000000e+13 ; 
gaus2.x0 = 8.840513e+15 ; gaus2.y0 = -1.015534e+08 ;  gaus2.sigma_x = -4.276958e+13 ; gaus2.sigma_y = 9.246286e+08 ;  gaus2.A = -2.429752e+19 ; 
gaus3.x0 = -4.241396e-01 ; gaus3.y0 = -6.386492e-07 ;  gaus3.sigma_x = 4.743911e-01 ; gaus3.sigma_y = 2.069032e-01 ;  gaus3.A = 1.123971e+30 ; 
gaus4.x0 = -4.217083e-01 ; gaus4.y0 = 3.842409e-07 ;  gaus4.sigma_x = 2.718323e-01 ; gaus4.sigma_y = 1.447877e-01 ;  gaus4.A = 3.725053e+30 ; 

 //// For I_Pulse2_30*Rho_e_30   TIME::
gaus1.x0 = 1.000000e-01 ; gaus1.y0 = -0.000000e+00 ;  gaus1.sigma_x = 2.400000e+00 ; gaus1.sigma_y = 4.000000e-01 ;  gaus1.A = 2.000000e+13 ; 
gaus2.x0 = -2.990009e+07 ; gaus2.y0 = 3.197690e+06 ;  gaus2.sigma_x = -1.578683e+08 ; gaus2.sigma_y = -6.780406e+06 ;  gaus2.A = 1.866019e+26 ; 
gaus3.x0 = -4.214189e-01 ; gaus3.y0 = -3.365750e-07 ;  gaus3.sigma_x = 4.656026e-01 ; gaus3.sigma_y = 2.035654e-01 ;  gaus3.A = 2.924066e+29 ; 
gaus4.x0 = -4.399315e-01 ; gaus4.y0 = 2.752652e-07 ;  gaus4.sigma_x = -2.703740e-01 ; gaus4.sigma_y = 1.524959e-01 ;  gaus4.A = 9.871485e+29 ; 

 //// For I_Pulse2_35*Rho_e_35   TIME::
gaus1.x0 = -1.500000e+00 ; gaus1.y0 = -0.000000e+00 ;  gaus1.sigma_x = 2.400000e+00 ; gaus1.sigma_y = 4.000000e-01 ;  gaus1.A = 2.000000e+10 ; 
gaus2.x0 = 3.480520e+04 ; gaus2.y0 = -1.060810e+06 ;  gaus2.sigma_x = 1.570583e+06 ; gaus2.sigma_y = 8.167302e+05 ;  gaus2.A = 3.882467e+25 ; 
gaus3.x0 = -4.284291e-01 ; gaus3.y0 = 7.610975e-07 ;  gaus3.sigma_x = -2.700529e-01 ; gaus3.sigma_y = 1.486990e-01 ;  gaus3.A = 1.039513e+29 ; 
gaus4.x0 = -4.110616e-01 ; gaus4.y0 = -3.320112e-06 ;  gaus4.sigma_x = 4.690216e-01 ; gaus4.sigma_y = 2.080050e-01 ;  gaus4.A = 2.638904e+28 ; 


'''

# guess_prms = [(-2.8, 1.5, 1, 0.5, 2.0e13),
#               (-2.8, -0.5, 1, 0.5, 2.0e13),
#               (-1.7, 0.5, 0.4, 1.5, 1.0e12),
#               (1.7, 2.0, 0.4, 1.5, 1.0e12),
#               (-2.7, -2.0, 0.4, 1.5, 1.0e13)
#              ]

# guess_prml = [(-2.0e+2,  9.1e+1, -2.0e+1, -1.29e+1
# , -4.39169913e+10, -2.33478842e+10)]

# ,(-2.66150850e+20,  9.49579891e+10, -2.06357311e+10, -1.88798129e+10
# , -4.39169913e+10, -2.33478842e+01)] 
                                     


'''
 -2.66150850e+21  9.49579891e+15 -2.06357311e+14 -1.88798129e+11
 -4.39169913e+11 -2.33478842e+07

'''

## How many functions were fitted
# gaus_no = len(guess_prms)
# loren_no = len(guess_prml)






## FITTED PARAMS FOR GAUSSIAN ENERGY DEPOSIION
'''
Fitted parameters:
[-5.09472709e-01  4.36820539e-01  5.01536752e-01  2.54298856e-01  1.99007486e+13 
-5.09467920e-01 -4.36822323e-01  5.01532148e-01  2.54297070e-01  1.99007062e+13  
1.24061914e+00 -1.62765076e-05  6.22783729e-01 -8.45522591e-01  6.58660628e+12 
-2.05679838e+00 -4.66312305e-05  4.66166971e-01  8.83421780e-01  6.10485144e+12]

'''

##---------------
## INITIAL GUESS DONE
##----------------------------------------------------------------------------

# We need to ravel the meshgrids of X, Y points to a pair of 1-D arrays.
xdata = np.vstack((X.ravel(), Y.ravel()))
# Do the fit, using our custom _gaussian function which understands our
# flattened (ravelled) ordering of the data points.
popt, pcov = curve_fit(_fun_fit, xdata, Z.ravel(), allpars , maxfev=600000)




# # The two-dimensional domain of the fit.
xmin, xmax, nx = -2.6, 2.6, 150
ymin, ymax, ny = -9.8, 9.8, 150
x, y = np.linspace(xmin, xmax, nx), np.linspace(ymin, ymax, ny)
X, Y = np.meshgrid(x, y)



fit = np.zeros(X.shape)
for i in range(gaus_no):
    fit += gaussian(X, Y, *popt[i*5:i*5+5])

for m in range(loren_no):
    # m = i + gaus_no 
    # print("inside loren")
    fit += lorentzian(X, Y, *popt[(gaus_no*5+m*6):((gaus_no*5+m*6)+6)])
    
for n in range(expo_no): # // gives quotient
   fit = fit * exponent(X, Y, *popt[(loren_no*6+gaus_no*5+n*4):((loren_no*6+gaus_no*5+n*4)+4)]) # pointer used

print('Fitted parameters:')
print(popt)



## The output required 
'''
gaus1.A = 1.99007486e+13 ; gaus1.x0 = -5.09472709e-01 ; gaus1.y0 = 4.36820539e-01; gaus1.sigma_x = 7.09271093e-01; gaus1.sigma_y = 3.59628433e-01 ;
gaus2.A = 1.99007062e+13 ; gaus2.x0 = -5.09467920e-01 ; gaus2.y0 =-4.36822323e-01; gaus2.sigma_x = 7.09275248e-01 ; gaus2.sigma_y =  3.59631310e-01 ;
gaus3.A =6.58660628e+12 ; gaus3.x0 = 1.24061914e+00; gaus3.y0 = 2.39198134e-05 ; gaus3.sigma_x = 8.80728819e-01; gaus3.sigma_y = 1.19575633e+00 ;
gaus4.A =6.10485144e+12 ; gaus4.x0 = -2.05679838e+00; gaus4.y0 =7.13013373e-06 ; gaus4.sigma_x = 6.59353934e-01; gaus4.sigma_y = 1.24931801e+00;


'''
dtime = datetime.datetime.now()
dtime = str(dtime)
outfile = open(fout_for_c, 'w+')
print("\n //// For "+ fname+'*'+fname2+"\t TIME::")
outfile.write("\n //// For "+ fname+'*'+fname2+"\t TIME::")
outfile.write(dtime)
for i in range(gaus_no):
    k = i+1
    paras =  popt[i*5:i*5+5] 
    print("gaus%d.x0 = %e ;" %(k, paras[0]), 'gaus%d.y0 = %e ; ' % (k, paras[1]),
          "gaus%d.sigma_x = %e ;" %(k, paras[2]), 'gaus%d.sigma_y = %e ; ' % (k, paras[3]), 'gaus%d.A = %e ; ' % (k, paras[4]) )

    outfile.write("\n gaus%d.x0 = %e ; \t gaus%d.y0 = %e ;  \t gaus%d.sigma_x = %e ; \t gaus%d.sigma_y = %e ;  \t gaus%d.A = %e ; "
     %(k, paras[0],k, paras[1], k, paras[2],k, paras[3],k, paras[4]) )

#  %(k, paras[0])
#  % (k, paras[1])
#  %(k, paras[2])
# %(k, paras[3])

    # outfile.write("gaus%d.x0 = %f ; " %i % paras[0]"gaus%d.y0 = %f ; " %i %paras[1])
    # np.savetxt(fout_for_c,"gaus1.x0 = %lf",paras[0] ) 

outfile.close() #Close the file when we’re done!    

# x,y,z equal sized 1D arrays


# Plot the 3D figure of the fitted function and the residuals.
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(X, Y, fit, cmap='plasma')
# cset = ax.contourf(X, Y, Z-fit, zdir='z', offset=-4.0e12, cmap='plasma')
ax.set_title("From Function")
ax.set_zlim(0,np.max(Z)+2)
plt.show()


# Plot the test data as a 2D image and the fit as overlaid contours.
fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow(Z, origin='bottom', cmap='plasma',
          extent=(xmin, xmax, ymin, ymax))
ax.contour(X, Y, fit,7, colors='b')
plt.show()