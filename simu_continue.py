'''
Created on 4 dec. 2017

@author: remi
'''
import sys
from qutip import *
import numpy as np
import pickle
from collections import namedtuple
from math import sqrt
import os
from multiprocessing import Process, Manager
import multiprocessing
import Averaging_discret
import time
Simu_continue = namedtuple("Simu_continue", ["epsilon", "alpha"])
pi=np.pi
up = basis(2, 0)
down=basis(2, 1)
def compute_continue(epsilon,alpha,times,E,psi0=basis(2, 0)):
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

    #print(H1_coeffR(0.96/epsilon,None))
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
    res=compute_continue(epsilon,alpha,times,E)
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

def get_times_d(E,alpha,epsilon):
    f_freq=Averaging_discret.Averaging_discret_sphere.get_f_freq(E,alpha)
    averag=Averaging_discret.Averaging_discret_sphere(f_freq=f_freq)
    averag.get_timer(epsilon)
    return averag.times

def prepare_continue_mpi(alpha,E,dt,epsilon,i):
    times_d=get_times_d(E,alpha,epsilon)
    ti=int(times_d[i]/(dt*epsilon))
    tf=int(times_d[i+1]/(dt*epsilon))
    times=np.arange(ti,tf,1)*dt
    #print(times)
    ez=up
    ex=(1/sqrt(2)*down+1/sqrt(2)*up)
    resultz=compute_continue(epsilon,alpha,times,E,psi0=ez)
    resultx=compute_continue(epsilon,alpha,times,E,psi0=ex)
    rx=np.array([resultx[0][-1],resultx[1][-1],resultx[2][-1]])
    rx=rx/np.linalg.norm(rx)
    rz=np.array([resultz[0][-1],resultz[1][-1],resultz[2][-1]])
    rz=rz/np.linalg.norm(rz)
    ry= np.cross(rz,rx)
    return np.transpose(np.array([rx,ry,rz]))

def prepare_continue_single_thread(alpha,E,dt,epsilon,path,i):
    tf=int(2*pi/epsilon)
    times=np.arange(0,tf,dt)
    res=compute_continue(epsilon,alpha,times,E)
    time.sleep(50*np.random.random(1))#eviter que tout le monde sauvegarde en meme temps
    times_d=get_time_d(E,alpha,epsilon)
    ## we load the continue simulation
    [x_c,y_c,z_c]=res
    ## we set the same time for both
    [x_c_p,y_c_p,z_c_p]=[np.interp(times_d, times*epsilon, x_c),np.interp(times_d, times*epsilon, y_c),np.interp(times_d, times*epsilon, z_c)]
    #we save the new result
    res=[x_c_p,y_c_p,z_c_p]
    save(path,res)

def save(path,res):
    try:
        with open(path, 'wb') as fp:
            pickle.dump(res, fp)
        os.system('scp '+path+ ' remi@129.104.219.230:test/')
        print("done")
    except:
        time.sleep(1)
        os.system('echo "{:d}" >> res/failed.txt'.format(i))
        os.system('rm '+path) 

#compute_continue(0.01,0,None,2)
# epsilon=0.1
# E=2
# alpha=0
# dt=0.001
# lk=len(get_times_d(E,alpha,epsilon))
# flot=np.eye(3)
# times_d=get_times_d(E,alpha,epsilon)
# print(times_d)
# t=int(times_d[lk-1]/(dt*epsilon))
# for k in range (lk-1):
#     flot=np.dot(prepare_continue_mpi(alpha,E,dt,epsilon,k),flot)
# print(flot)
# r=np.dot(flot,[0,0,1])
# r=r/np.linalg.norm(r)
# print(r)
# tf=int(2*pi/(epsilon*dt))*dt
# times=np.arange(0,tf,dt)
# res=compute_continue(epsilon,alpha,times,E)
# res=np.array(res)
# print(res[:,t])
# print(res[:,-1])
# print(len(res[0]))