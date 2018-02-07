import pickle
import Periodic_system
import Averaging_discret
import matplotlib.pyplot as plt
import numpy as np

with open('simu_rapide_xyztheta_500pt_-3.5dt', 'rb') as fp:
    xyztheta = pickle.load(fp)

[lst_x,lst_y,lst_z,lst_theta]=xyztheta
x=lst_x[2]
y=lst_y[2]
z=lst_z[2]
theta=lst_theta[2]

perio=Periodic_system.Periodic_system(False)
perio.regularize(theta)
#perio.plot_result(x,y,z,theta)

f_freq=Averaging_discret.Averaging_discret_sphere.get_f_freq(2,-0.5)
f_axis,f_theta=Averaging_discret.Averaging_discret_sphere.get_f_axis_theta(x,y,z,theta)
#times=np.linspace(0,2*np.pi,989)
#fig=plt.figure()
#plt.plot(times,f_freq(times))
#plt.plot(times,f_axis(times)[0])
#plt.plot(times,f_axis(times)[1])
#plt.plot(times,f_axis(times)[2])
#plt.plot(times,f_theta(times))
#plt.show()
averag=Averaging_discret.Averaging_discret_sphere(f_theta,f_axis,f_freq)
#print(len(averag.get_z(0.001)))
averag.plot_multi_N([0.1,0.01,0.001])
#times,xyz,theta=averag.get_z(0.01)#

#print(times)
#time=np.linspace(0,2*np.pi,len(times))
#fig=plt.figure()
#plt.plot(time,times)
#plt.plot(times,xyz[:,0])
#plt.plot(times,xyz[:,1])
#plt.plot(times,xyz[:,2])
#plt.plot(times,f_theta(times))
#plt.show()