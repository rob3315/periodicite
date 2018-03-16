import unittest
from simu_continue import *
import numpy as np

class test_simu_continue(unittest.TestCase):
	def setUp(self):
		self.E=2
		self.alpha=0

	def test_mpi_result(self):
		epsilon=0.1
		dt=0.001
		lk=len(get_times_d(self.E,self.alpha,epsilon))
		flot=np.eye(3)
		times_d=get_times_d(self.E,self.alpha,epsilon)
		t=int(times_d[lk-1]/dt)
		for k in range (lk-1):
			flot=np.dot(prepare_continue_mpi(self.alpha,self.E,dt,epsilon,k),flot)
		res_mpi=np.dot(flot,[0,0,1])
		tf=int(2*pi/epsilon)
		times=np.arange(0,tf,dt)
		res=compute_continue(epsilon,self.alpha,times,self.E)
		print(res)
		res_continue=np.array([res[0][t],res[1][t],res[2][t]])
		print(res_continue)
		print('t',t,tf)
		print(flot)
		self.assertAlmostEqual(np.linalg.norm(res_continue-res_mpi),0,places=4)


if __name__=='__main__':
	unittest.main()