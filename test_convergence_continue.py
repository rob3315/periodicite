'''
Created on 4 dec. 2017

@author: remi
'''
from qutip import *
import numpy as np
import pickle

pi=np.pi
up = basis(2, 0)
down=basis(2, 1)

def compute(epsilon,alpha,times):
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
    HC=[H0,[sigmax(),H1_coeffCx],[sigmay(),H1_coeffCy]]

    resultR = mesolve(HR, psi0, times, [], [sigmaz()])
    resultC = mesolve(HC, psi0, times, [], [sigmaz()])
    #print(resultR.expect[0][-1])
    return(resultR.expect[0][-1],resultC.expect[0][-1])


lst_alpha=np.arange(-1,1,0.1)
E=2
dt=0.001
#dt=10
lst_epsilon=np.array([0.05,0.02,0.01,0.008,0.005,0.002,0.001])
resR=np.zeros((len(lst_epsilon),len(lst_alpha)))
resC=np.zeros((len(lst_epsilon),len(lst_alpha)))
i=0
for epsilon in lst_epsilon:
    tf=int(2*pi/epsilon)
    times=np.arange(0,tf,dt)
    #compute(epsilon,0,times)
    def func1(alpha): return compute(epsilon,alpha,times)
    a,b = parfor(func1, lst_alpha)
    print(a,b)
    resR[i,:]=a
    resC[i,:]=b
    i+=1

with open('outfileR', 'wb') as fp:
    pickle.dump(resR, fp)
with open('outfileC', 'wb') as fp:
    pickle.dump(resC, fp)
print("done")
