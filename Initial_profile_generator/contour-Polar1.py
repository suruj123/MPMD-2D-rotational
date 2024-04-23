#! /usr/bin/env python
from pylab import *
from numpy import *
import matplotlib.patches as patches
from matplotlib import pyplot as plt
import matplotlib.patches as patches

import matplotlib as mpl
mpl.rc('figure', max_open_warning = 0)
a = 1
b = 45*45
z = zeros(b)
VX = zeros(b)
VY = zeros(b)
#data1 = loadtxt("STATE1")
for i in range(0,1,1):
   data = loadtxt('QUIVER_time_1_'+str(i))
   #data = loadtxt('QUIVER_time_1000')
   x =  data[:,0]
   y =  data[:,1]
   vx = data[:,2]
   vy = data[:,3]
   ke =  data[:,5]
   y1 = y
   #VX += vx
   #VY += vy
   #z += Gamma

#Gamma = z/400
   T = ke
   xi = np.linspace(min(x), max(x))
   yi = np.linspace(min(y), max(y))
   X, Y = np.meshgrid(xi, yi)
   T1 = griddata(x, y, T, xi, yi,'linear')
   #values = values.reshape(len(zeniths), len(azimuths))
   #ke = np.array(ke)
   #ke = ke.reshape(len(x), len(y))

  
   
 
   fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
   #ax.set_rlim(50,110)
   ax.set_rorigin(-20)

   #contourplot = ax.contourf(theta, r, values)
   contourplot = ax.contourf(Y, X, ke)
   #ax.set_yticklabels([])
   plt.colorbar(contourplot, shrink=1.0, pad=0.08)





   #cbar.ax.tick_params( labelsize=10)#,color="black", width=1, ,)
   #cbar.set_label(r"$ T = 1/\Gamma$", fontsize=22)
   #ax.set_aspect('1.0')
  
   plt.title( r'$\Gamma = 100$, time = %lf T' % (float(i)/170),fontsize=10, pad= 12.0)
   
   #xlabel(r"$x$", fontsize=10)
   #ylabel(r"$y$",fontsize=10)
   #ax.tick_params(axis='both', which='major', labelsize=22)
   #ax.yaxis.set_major_locator(plt.MultipleLocator(40))
   #ax.xaxis.set_major_locator(plt.MultipleLocator(60))


   #quiver(x, y, vx, vy,color='white')#,scale=0.08, units='inches',width = 0.025)
   savefig("Vorticitya"+str(i)+".png")
   #show()
   #savefig("AR2_Contour_Temp_RTD_1p4.pdf")
   #plt.show()   
   
#show()

