from simu_continue import *
lst_alpha=np.arange(-1,1,0.1)
path="res/simu_continue_compact_{:f}_{:f}_-4dt"
E=2
dt=0.0001
lst_epsilon=np.array([0.05,0.02,0.01,0.008,0.005,0.002,0.001,0.0005,0.0001])
i=int((sys.argv[1]))
print(i)
epsilon=lst_epsilon[i//len(lst_alpha)]
alpha=lst_alpha[i % len(lst_alpha)]
path=path.format(epsilon,alpha)
prepare_continue_single_thread(alpha,E,dt,epsilon,path,i)
#prepare_continue(lst_alpha,E,dt,lst_epsilon,path)