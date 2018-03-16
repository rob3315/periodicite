from qutip import *
import matplotlib.pyplot as plt
import numpy as np
#from mpl_toolkits.mplot3d import Axes3D
pi=np.pi
from math import sqrt
import pickle
import Periodic_system
import Averaging_discret

lst_alpha=np.arange(-1,1,0.1)
#lst_epsilon=np.array([0.05,0.02,0.01,0.008,0.005])#,0.002,0.001])
lst_epsilon=[0.1,0.05,0.01,0.005,0.001,0.0005,0.0001,0.00005]
#epsilon=lst_epsilon[5]
E=2
i=10
alpha=int(lst_alpha[i]*10)/10
#dt=0.001

# we load and compute the discret simulation
#with open('simu_xyztheta_1000pt_-4dt', 'rb') as fp:
#    xyztheta = pickle.load(fp)
#[lst_x_s,lst_y_s,lst_z_s,lst_theta]=xyztheta#
#x_s=lst_x_s[i]
#y_s=lst_y_s[i]
#z_s=lst_z_s[i]
#theta=lst_theta[i]
#print('theta')
with open('res/simu_xyztheta_1000pt-5dt_-{:f}'.format(alpha), 'rb') as fp:
    xyztheta = pickle.load(fp)
[x_s,y_s,z_s,theta]=xyztheta##
x_s=np.array(x_s)
x_s=x_s.flatten()
y_s=np.array(y_s)
y_s=y_s.flatten()
z_s=np.array(z_s)
z_s=z_s.flatten()
theta=np.array(theta)
theta=theta.flatten()


#perio=Periodic_system.Periodic_system(False)
#perio.regularize(theta)
#perio.plot_result(x,y,z,theta)

f_freq=Averaging_discret.Averaging_discret_sphere.get_f_freq(2,alpha)
f_axis,f_theta=Averaging_discret.Averaging_discret_sphere.get_f_axis_theta(x_s,y_s,z_s,theta)
averag=Averaging_discret.Averaging_discret_sphere(f_theta,f_axis,f_freq)

lx=len(lst_epsilon)//2 +len(lst_epsilon)%2
fig1, axs = plt.subplots(lx, 2, sharex=True, sharey=True)

for k in range(len(lst_epsilon)):
	epsilon=lst_epsilon[k]
	ax=axs[k%lx,k//lx]
	[xp,yp,zp]=np.transpose(averag.get_z(epsilon))
	times_d=averag.times
	## we load the continue simulation
	path="res/simu_continue_compact_{:f}_{:f}_-5dt"
	path=path.format(epsilon,alpha)
	print(path)
	with open(path,'rb') as fp:
	    res_c=pickle.load(fp)
	[x_c,y_c,z_c]=res_c
	tf=int(2*pi/epsilon)
	#times_c=np.arange(0,tf,dt)
	#print(len(x_c),len(times_c))

	## we set the same time for both
	#[x_c_p,y_c_p,z_c_p]=[np.interp(times_d, times_c*epsilon, x_c),np.interp(times_d, times_c*epsilon, y_c),np.interp(times_d, times_c*epsilon, z_c)]

	[x_c_p,y_c_p,z_c_p]=res_c
	#plot

	#plt.plot(times_c*epsilon,y_c)
	ax.plot(times_d,zp,'gold')#
	ax.plot(times_d,z_c_p,'black')#
	ax.plot(times_d,yp,'r')#
	ax.plot(times_d,y_c_p,'y')#
	ax.plot(times_d,xp,'b')#
	ax.plot(times_d,x_c_p,'g')#

	ax.set_title('epsilon = {:f}'.format(epsilon))
print(axs)
axs[0][0].legend(['z discret','z continue','y discret','y continue','x discret','x continue'])

fig2 = plt.figure()
ax2 = plt.gca()
averag.plot_final(epsilon,ax2)
plt.show()
#print(times_c)
#print(times_d)