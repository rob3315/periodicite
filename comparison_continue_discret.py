from qutip import *
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
pi=np.pi
from math import sqrt
import pickle
import Periodic_system
import Averaging_discret

lst_alpha=np.arange(-1,1,0.1)
lst_epsilon=np.array([0.05,0.02,0.01,0.008,0.005,0.002,0.001])
epsilon=lst_epsilon[5]
E=2
i=10
alpha=lst_alpha[i]
dt=0.001

## we load and compute the discret simulation
with open('simu_xyztheta_1000pt_-4dt', 'rb') as fp:
    xyztheta = pickle.load(fp)

[lst_x_s,lst_y_s,lst_z_s,lst_theta]=xyztheta

x_s=lst_x_s[i]
y_s=lst_y_s[i]
z_s=lst_z_s[i]
theta=lst_theta[i]

perio=Periodic_system.Periodic_system(False)
perio.regularize(theta)
#perio.plot_result(x,y,z,theta)

f_freq=Averaging_discret.Averaging_discret_sphere.get_f_freq(2,alpha)
f_axis,f_theta=Averaging_discret.Averaging_discret_sphere.get_f_axis_theta(x_s,y_s,z_s,theta)
averag=Averaging_discret.Averaging_discret_sphere(f_theta,f_axis,f_freq)
[xp,yp,zp]=np.transpose(averag.get_z(epsilon))
times_d=averag.times
## we load the continue simulation
path="res_continue/simu_continue_{:f}_{:f}_-4dt"
path=path.format(epsilon,alpha)
with open(path,'rb') as fp:
    res_c=pickle.load(fp)
[x_c,y_c,z_c]=res_c
tf=int(2*pi/epsilon)
times_c=np.arange(0,tf,dt)
print(len(x_c),len(times_c))

## we set the same time for both
[x_c_p,y_c_p,z_c_p]=[np.interp(times_d, times_c*epsilon, x_c),np.interp(times_d, times_c*epsilon, y_c),np.interp(times_d, times_c*epsilon, z_c)]
print(len(times_d))
#plot
fig=plt.figure()
#plt.plot(times_c*epsilon,y_c)
plt.plot(times_d,yp)#
plt.plot(times_d,y_c_p)#
plt.plot(times_d,xp)#
plt.plot(times_d,x_c_p)#
plt.plot(times_d,zp)#
plt.plot(times_d,z_c_p)#

##plt.plot(times_d,y_c_p-yp)
##plt.plot(times_d,z_c_p-zp)
plt.show()
#print(times_c)
#print(times_d)