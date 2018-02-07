'''
Created on 4 dec. 2017

@author: remi
'''
import sys
from qutip import *
import numpy as np
import pickle
from collections import namedtuple
import os
from multiprocessing import Process, Manager
import multiprocessing
import time
Simu_continue = namedtuple("Simu_continue", ["epsilon", "alpha"])
pi=np.pi
up = basis(2, 0)
down=basis(2, 1)
def compute_continue(epsilon,alpha,times):
    psi0=basis(2, 0)
    H0=(E+alpha)*sigmaz()

    u1 = lambda t : 1.5*(1-np.cos(t))
    u2 = lambda t : -3 *np.cos(t/2)
    U2= lambda t: -6 * np.sin(t/2) # primitive of u2 such as U2(0)=0

    def H1_coeffR(t,args):
        return u1(epsilon*t) *1/2 *np.real( np.exp( -1j*(2*E*t +  U2(epsilon*t)/epsilon ) ))

    def H1_coeffCx(t,args):
        return u1(epsilon*t) *1/4*np.real( np.exp( -1j*(2*E*t +  U2(epsilon*t)/epsilon ) ))

    def H1_coeffCy(t,args):
        return -u1(epsilon*t) *1/4*np.imag( np.exp( -1j*(2*E*t +  U2(epsilon*t)/epsilon ) ))


    HR=[H0,[sigmax(),H1_coeffR]]
    #HC=[H0,[sigmax(),H1_coeffCx],[sigmay(),H1_coeffCy]]

    resultR = mesolve(HR, psi0, times, [], [sigmax(),sigmay(),sigmaz()])
    #resultC = mesolve(HC, psi0, times, [], [sigmaz()])
    #print(resultR.expect[0][-1])
    #return([epsilon,alpha,times])
    return(resultR.expect)


def calculate(func, args):
    print('eps,alpha,dic',args)
    result = func(*args)
    return result
    #return '%s did epsilon = %f alpha = %f ' % (multiprocessing.current_process().name,func.__name__, args            )
def calculatestar(args):
    #print('hi2',args)
    return calculate(args[0],args[1])

def f_aux(epsilon,alpha,dic):
    tf=int(2*pi/epsilon)
    times=np.arange(0,tf,dt)
    res=compute_continue(epsilon,alpha,times)
    simu=Simu_continue(epsilon,alpha)
    #print(dic)
    dic[simu]=res
    return None
    #return dic

def prepare_continue_multi_thread(lst_alpha,E,dt,lst_epsilon,path):
    manager=Manager()
    dicR=manager.dict()



    TASKS=[k for i in [[(f_aux,(epsilon,alpha,dicR)) for epsilon in lst_epsilon] for alpha in lst_alpha] for k in i]
    #print(TASKS)
    PROCESSES=4
    pool = multiprocessing.Pool(PROCESSES)
    imap_it = pool.map(calculatestar, TASKS)
    print('imap',imap_it)

    print(dicR)
    with open(path, 'wb') as fp:
        pickle.dump(dicR, fp)
    print("done")

def prepare_continue_single_thread(alpha,E,dt,epsilon,path,i):
    tf=int(2*pi/epsilon)
    times=np.arange(0,tf,dt)
    res=compute_continue(epsilon,alpha,times)
    time.sleep(50*np.random.random(1))#eviter que tout le monde saauvegarde en meme temps
    try:
        with open(path, 'wb') as fp:
            pickle.dump(res, fp)
        os.system('scp '+path+ ' remi@129.104.219.230:test/')
        print("done")
    except:
        time.sleep(1)
        os.system('echo "{:d}" >> res/failed.txt'.format(i))
        os.system('rm '+path) 

lst_alpha=np.arange(-1,1,0.1)
path="res/simu_continue_{:f}_{:f}_-4dt"
E=2
dt=0.001
lst_epsilon=np.array([0.05,0.02,0.01,0.008,0.005,0.002,0.001])
i=int((sys.argv[1]))
print(i)
epsilon=lst_epsilon[i//len(lst_alpha)]
alpha=lst_alpha[i % len(lst_alpha)]
path=path.format(epsilon,alpha)
prepare_continue_single_thread(alpha,E,dt,epsilon,path,i)
#prepare_continue(lst_alpha,E,dt,lst_epsilon,path)
