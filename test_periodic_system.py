import unittest
from Periodic_system import *
import numpy as np

class test_periodic_system(unittest.TestCase):
	def setUp(self):
		self.E=2
		self.alpha=0

	def test_get_axis_angle_aux(self):
		perSys=Periodic_system(False)
		for l in range(1000):
			k=2*np.random.rand(3)-np.ones(3)
			k=k/np.linalg.norm(k)
			theta=np.pi*np.random.rand(1)
			I=np.eye(3)
			M=np.zeros((3,3))
			for i in range(3):
				v=I[:,i]
				M[:,i]=v*np.cos(theta)+np.cross(k,v)*np.sin(theta)+np.dot(k,v)*(1-np.cos(theta))*k
			print(k,theta)
			print(M)
			kk,theta2=perSys.get_axis_angle_aux(M)
			print(kk,theta2)
			self.assertAlmostEqual(np.linalg.norm(kk-k),0,places=6)
			self.assertAlmostEqual(theta2,theta[0],places=4)



if __name__=='__main__':
	unittest.main()