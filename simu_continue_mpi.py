from mpi4py import MPI
from simu_continue import *
nproc_u=int(sys.argv[2])-1 #number of processeur processing
main=nproc_u # main proc number
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
epsilon=float((sys.argv[1]))
print(epsilon)
E=2
alpha=-0.5
dt=0.00001
path= "res_mpi/simu_continue_compact_{:f}_{:f}_-5dt".format(epsilon,alpha)
lk=len(get_times_d(E,alpha,epsilon))
if rank == main:
	#flot=prepare_continue_mpi(alpha,E,dt,epsilon,0)
	times_d=get_times_d(E,alpha,epsilon)
	lst_pos=np.zeros((3,lk))
	lst_pos[:,0]=np.array([0,0,1])
	#lst_pos[:,1]=np.dot(flot,lst_pos[:,0])
	for k in range (0,lk-1):
		#recvarray = np.zeros((3,3),dtype=np.float64)
		recvarray=comm.recv( source = k%nproc_u)
		lst_pos[:,k+1]=np.dot(recvarray,lst_pos[:,k])
		lst_pos[:,k+1]=lst_pos[:,k+1]/np.linalg.norm(lst_pos[:,k+1])
	print(lst_pos)
	print(len(lst_pos[0]))
	save(path,lst_pos)
else:
	i=rank
	while i<lk:
		data=prepare_continue_mpi(alpha,E,dt,epsilon,i)
		#print(data.dtype)
		comm.send(data, dest = main)
		i=i+nproc_u
	print("{:d} done".format(rank))