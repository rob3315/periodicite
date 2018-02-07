import numpy as np
import matplotlib.pyplot as plt
from abc import ABC, abstractmethod
class Averaging_discret(ABC):
	def __init__(self,f_theta):
		self.f_theta=f_theta
		super().__init__()

	@abstractmethod
	def rotate(self,angle,axis,z):
		pass

	@abstractmethod
	def single_plot_z(self,n,z,ax):
		pass

	@abstractmethod
	def get_z(self,n):
		pass

	#@abstractmethod
	def plot_final(self,n):
		pass

	def compute(self,lst_theta,lst_axis,z0,z):
		"""compute the discret averaging from z0 and store it in z
		rotate must be implemented"""
		z[0]=z0
		for k in range(len(lst_theta)-1):
			z[k+1]=self.rotate(lst_theta[k],lst_axis[k],z[k])


	def plot_multi_N(self,lst_n):
		fig, axs = plt.subplots(len(lst_n)+1, 1, sharex=True, sharey=False)

		for k in range(len(lst_n)):
			n=lst_n[k]
			z=self.get_z(n)
			self.single_plot_z(n,z,axs[k])
		self.plot_final(n,axs[-1])
		plt.show()

class Averaging_discret_flat(Averaging_discret):
	def __init__(self,f_theta):
		super().__init__(f_theta)

	def rotate(self,angle, axis,z):
		return np.exp(1j*angle)*(z-axis)+axis

	def get_z(self,n):
		times=np.arange(n)
		lst_theta=self.f_theta(times/n)
		lst_axis=np.array(times/n,dtype=complex)
		z=np.zeros(len(lst_theta),dtype=complex)
		self.compute(lst_theta,lst_axis,0+0j,z)
		return z
	def single_plot_z(self,n,z,ax):
		ax.plot(np.real(z), np.imag(z))
		ax.legend(["n = "+str(n)])

	def plot_final(self,n,ax):
		times=np.arange(n)/n
		lst_theta=np.real(self.f_theta(times))
		ax.plot(times,lst_theta)
		ax.plot(times,np.zeros(len(times)),'orange')
		ax.legend(["angle","zeros"])

class Averaging_discret_sphere(Averaging_discret):
	def get_f_freq(E,alpha):
		#control_u1 = lambda t : 1.5*(1-np.cos(t))
		control_u2 = lambda t : -3 *np.cos(t/2)
		return lambda t : (2*E+control_u2(t))/(2*np.pi)
	def get_f_axis_theta(x,y,z,theta):
		n=len(x)
		times=np.linspace(0,2*np.pi,n)
		return (lambda t : [np.interp(t, times, x),np.interp(t, times, y),np.interp(t, times, z)], lambda t : np.interp(t, times, theta))

	def __init__(self,f_theta,f_axis=None,f_freq=None):
		super().__init__(f_theta)
		self.f_axis=f_axis
		self.f_freq=f_freq

	def rotate(self,angle, axis,z):
		"""Rodrigues' rotation formula
		axis must be normalized"""
		#print(z,axis,angle)
		axis=axis/np.linalg.norm(axis) # we normalize anyway because of numerical error
		return np.cos(angle)*z + np.sin(angle) *np.cross(axis,z)+ (1-np.cos(angle))*(np.dot(axis,z))*axis
	def get_timer(self,epsilon):
		t=0
		self.times=[t] # between 0 and 2*pi
		while t<int(2*np.pi/epsilon):
			deltat=1/self.f_freq(t*epsilon)
			deltat=(1/self.f_freq(t*epsilon)+1/self.f_freq((t+deltat)*epsilon))/2
			t=t+deltat # irregular paces
			self.times.append(epsilon*t) # normalize times
		self.times=np.array(self.times)

	def get_z(self,epsilon):
		self.get_timer(epsilon)
		self.lst_theta=self.f_theta(self.times)
		self.lst_axis=np.transpose(np.array(self.f_axis(self.times)))
		z=np.zeros((len(self.lst_theta),3))
		
		#return(self.times,lst_axis,self.lst_theta)

		self.compute(self.lst_theta,self.lst_axis,(0,0,1),z)
		return z

	def single_plot_z(self,n,z,ax):
		ax.plot(self.times, z[:,0])
		ax.plot(self.times, z[:,1])
		ax.plot(self.times, z[:,2])
		ax.plot(self.times, z[:,0]**2+z[:,1]**2+z[:,2]**2)

		ax.legend(["x","y","z",'norme'])

	def plot_final(self,n,ax):
		ax.plot(self.times,np.mod(self.lst_theta,2*np.pi))
		ax.plot(self.times,self.lst_axis[:,0])
		ax.plot(self.times,self.lst_axis[:,1])
		ax.plot(self.times,self.lst_axis[:,2])
		ax.plot(self.times,self.lst_axis[:,0]**2+self.lst_axis[:,1]**2+self.lst_axis[:,2]**2)
		ax.legend(["angle","axis x","axis y","axis z","norm axis"])


#lst_theta=lambda times : np.pi/2*np.cos(2*np.pi*times,dtype=complex)
##lst_theta=np.pi/2*np.cos((1/3)*np.pi*times,dtype=complex)#np.sin(times,dtype=complex)
#test=Averaging_discret_flat(lst_theta)
#lst_n=[5,10,100,1000,10000]#,100000,1000000]
#test.plot_multi_N(lst_n)

