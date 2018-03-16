import os

lst_epsilon=[0.1,0.05,0.01,0.005,0.001,0.0005,0.0001,0.00005]
##lst_epsilon=[0.00001,0.000005]
nproc=1000
alpha=-0.6
for epsilon in lst_epsilon:
	os.system("source activate qutip-env \n salloc -n {:d} mpiexec python simu_continue_mpi.py {:f} {:d} -v".format(nproc,epsilon,nproc))
#os.system("source activate qutip-env \n salloc -n {:d} mpiexec python simu_continue_mpi.py {:f} -v".format(nproc,alpha))