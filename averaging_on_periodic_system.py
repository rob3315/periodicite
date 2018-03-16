import pickle
import Periodic_system
import Averaging_discret
import matplotlib.pyplot as plt
import numpy as np

#with open('simu_xyztheta_1000pt_-4dt', 'rb') as fp:
#    xyztheta = pickle.load(fp)
#lst_alpha=np.arange(-1,1,0.1)#

#i=4
#alpha=lst_alpha[i]
#[lst_x,lst_y,lst_z,lst_theta]=xyztheta
#x=lst_x[i]
#y=lst_y[i]
#z=lst_z[i]
#theta=lst_theta[i]

alpha=0
with open('res/simu_xyztheta_10000pt-5dt_-0.000000', 'rb') as fp:
    xyztheta = pickle.load(fp)
[x_s,y_s,z_s,theta]=xyztheta##
x=np.array(x_s)
x=x.flatten()
y=np.array(y_s)
y=y.flatten()
z=np.array(z_s)
z=z.flatten()
theta=np.array(theta)
theta=theta.flatten()


perio=Periodic_system.Periodic_system(False)
perio.regularize(theta)
#perio.plot_result(x,y,z,theta)

f_freq=Averaging_discret.Averaging_discret_sphere.get_f_freq(2,alpha)
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
averag.plot_multi_N([1E-2,1E-3,1E-4])#,1E-5])
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