from Periodic_system import *
from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


nb_task=500
E=2
dt=0.0001
by_task=500

test=Periodic_system(False)
lst_alpha=np.arange(-1,1,0.1)
res_tot=[]
i=int((sys.argv[1]))
alpha=lst_alpha[i]
ti=0.95
tf=0.98
if rank == 0:
	#pts=np.linspace(0,2*pi,nb_task)
	pts=np.linspace(ti,tf,nb_task)
else:
   pts = None
   
pts = comm.scatter(pts, root=0)
pts=np.linspace(pts,pts+(tf-ti)/(nb_task-1),by_task)

(x,y,z,theta)=test.launch_simu(nb_point=nb_task,E=E,dt=dt,c=1,alpha=alpha,pts=pts)
#(x,y,z,theta)=test.launch_simu(E=E,dt=dt,c=1,alpha=alpha,pts=[pts])

lst_x=None
lst_x=comm.gather(x,root=0)
res_tot.append(lst_x)

lst_y=None
lst_y=comm.gather(y,root=0)
res_tot.append(lst_y)

lst_z=None
lst_z=comm.gather(z,root=0)
res_tot.append(lst_z)

lst_theta=None
lst_theta=comm.gather(theta,root=0)
res_tot.append(lst_theta)
i
if rank == 0:
	try:
	    path='res/simu_xyztheta_zeros_{:d}_euc_pt-4dt_{:f}'.format(nb_task*by_task,alpha)
	    with open(path, 'wb') as fp:
	        pickle.dump(res_tot, fp)
	    os.system('scp '+path+ ' remi@129.104.219.230:test/')
	    #os.system('rm '+path) 
	    print("done")
	except:
	    time.sleep(1)
	    os.system('echo "{:d}" >> res/failed3.txt'.format(i))