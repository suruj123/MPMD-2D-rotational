#! /usr/bin/env python
from pylab import *
from numpy import *
import matplotlib.patches as patches
from matplotlib import pyplot as plt
import matplotlib.patches as patches

import matplotlib as mpl
mpl.rc('figure', max_open_warning = 0)
a = 1
b = 7056
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


   fig = figure(figsize=(20,10))
   ax = fig.add_subplot(111)
   
   
   
   #cp=plt.contour(X,Y,T1,100, levels=np.linspace(-0.4,0.4), extend='both',cmap=mpl.cm.jet)#,np.linspace(-1,1,11))
   #cp = contourf(X,Y,T1,100, extend='both',cmap=mpl.cm.jet)
   #cp = imshow(T1,cmap = cmap, extent =[X.min(), X.max(), Y.min(), Y.max()])#, aspect = "2.2")#, vmax= 0.6, vmin = -0.6)#, norm = norm)# interpolation = "spline16", aspect = "auto")

   #levels = np.linspace(-0.6, 0.6, 30)
   norm = cm.colors.Normalize(vmax=T1.max(), vmin=T1.min())
   #cmap = cm.seismic
   cmap = mpl.cm.jet
   cp = imshow(T1,cmap = cmap, extent =[X.min(), X.max(), Y.min(), Y.max()], vmax= 0.158339, vmin = -0.138204, norm = norm, interpolation = "spline16", aspect = "auto")
   cbar = plt.colorbar(cp,cmap = cmap, orientation='vertical')#, ticks = [0, 0.026, 0.045])
   
   Drawing_colored_circle = plt.Circle(( 0.0 , 0.0 ), 50, color = 'w' )
   ax.add_artist( Drawing_colored_circle )
   cbar.ax.tick_params( labelsize=10)#,color="black", width=1, ,)
   #cbar.set_label(r"$ T = 1/\Gamma$", fontsize=22)
   #ax.set_aspect('1.0')
  
   #plt.title( r'circulation = 3,$\Gamma = 100$,time ='+str(i),fontsize=22)
   plt.title( r'time ='+str(i),fontsize=22)
  # plt.xlim(-500,500)
  # plt.ylim(-500,500)	
   xlabel(r"$x$", fontsize=10)
   ylabel(r"$y$",fontsize=10)
   ax.tick_params(axis='both', which='major', labelsize=15)
   ax.yaxis.set_major_locator(plt.MultipleLocator(40))
   ax.xaxis.set_major_locator(plt.MultipleLocator(40))


   #quiver(x, y, vx, vy,color='white')#,scale=0.08, units='inches',width = 0.025)
   savefig("Vorticity"+str(i)+".jpg")
   #show()
   savefig("Vorticity.pdf")
   
   
   show()

