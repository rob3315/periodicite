import unittest
from Averaging_discret import *
import numpy as np

class testAveraging(unittest.TestCase):
	def setUp(self):
		self.avg_sphere=Averaging_discret_sphere(None,None,None)
	def testRotate(self):
		[e1,e2,e3]=np.eye(3)
		#direct basis
		self.assertAlmostEqual(np.linalg.norm(self.avg_sphere.rotate(np.pi/2,e1,e2)-e3),0)
		self.assertAlmostEqual(np.linalg.norm(self.avg_sphere.rotate(np.pi/2,e2,e3)-e1),0)
		self.assertAlmostEqual(np.linalg.norm(self.avg_sphere.rotate(np.pi/2,e3,e1)-e2),0)

		self.assertAlmostEqual(np.linalg.norm(self.avg_sphere.rotate(np.pi/2,e2,e1)+e3),0)
		self.assertAlmostEqual(np.linalg.norm(self.avg_sphere.rotate(np.pi/2,e3,e2)+e1),0)
		self.assertAlmostEqual(np.linalg.norm(self.avg_sphere.rotate(np.pi/2,e1,e3)+e2),0)

		#pi rotation
		self.assertAlmostEqual(np.linalg.norm(self.avg_sphere.rotate(np.pi,e2,e1)+e1),0)
		self.assertAlmostEqual(np.linalg.norm(self.avg_sphere.rotate(np.pi,e3,e2)+e2),0)
		self.assertAlmostEqual(np.linalg.norm(self.avg_sphere.rotate(np.pi,e1,e3)+e3),0)

		#
		for i in range(2):
			z1=np.random.random(3)
			z1=z1/np.linalg.norm(z1)
			for j in range(10000):
				z2=np.random.random(3)
				z2=z2/np.linalg.norm(z2)
				angle=3*np.pi*np.random.random(1)
				z1=self.avg_sphere.rotate(angle,z2,z1)
			self.assertAlmostEqual(np.linalg.norm(z1),1)
			print(np.linalg.norm(np.linalg.norm(z1)))


if __name__=='__main__':
	unittest.main()