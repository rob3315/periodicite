from qutip import *
import matplotlib.pyplot as plt
import numpy as np
#from mpl_toolkits.mplot3d import Axes3D
pi=np.pi
import sys
import os
import time
from math import sqrt
import pickle
class Periodic_system():
    qutip=False
    def __init__(self,plot):
        self.up = basis(2, 0)
        self.down=basis(2, 1)
        self.control_u1 = lambda t : 1.5*(1-np.cos(t))
        self.control_u2 = lambda t : -3 *np.cos(t/2)
        self.control_U2= lambda t: -6 * np.sin(t/2) # primitive of u2 such as U2(0)=0
        self.plot_res=plot # plot of the result
        self.plot_inter=False # plot of every configuration
    def compute(self,w,psi0,c,dt):
        times=np.arange(0,2*np.pi/w,dt)
        def H1_coeffR(t,args):
            return c*1/2*np.real( np.exp( -1j*(w*t) ))

        if Periodic_system.qutip:
            H0=(self.E+self.alpha)*sigmaz()
            HR=[H0,[sigmax(),H1_coeffR]]
            resultR = mesolve(HR, psi0, times, [], [sigmax(), sigmay(),sigmaz()])
            return(resultR.expect,H1_coeffR(times,()))
        else :
            #euler schema
            #psi0 = np.array([0,0,1])# TODO implement the conversion
            result=np.zeros((3,2))
            result[:,0]=psi0
            H=lambda t: np.array([[0,-2*(self.E+self.alpha),0],
                [2*(self.E+self.alpha),0,-2*H1_coeffR(t,())],
                [0,2*H1_coeffR(t,()),0]])
            for k in range(1,len(times)):
                result[:,1]=result[:,0]+np.dot(H(times[k]),result[:,0])*dt
                result[:,1]=result[:,1]/(np.linalg.norm(result[:,1]))
                result[:,0]=result[:,1]
            return(result,H1_coeffR(times,()))


    def get_axis_angle(self,dt,w,c):
    # we compute the evolution of ex, ey and ez
    # in fact it is enough to get the evolution of ex and ey because
    #the evolution matrix in in SO_3 
        if Periodic_system.qutip:
            psi0=(1/sqrt(2)*self.down+1/sqrt(2)*self.up)
        else :
            psi0 = np.array([1,0,0])
        resultR,control=self.compute(w,psi0,c,dt)
        rx=[resultR[0][-1],resultR[1][-1],resultR[2][-1]]

        if Periodic_system.qutip:
            psi0=(1j/sqrt(2)*self.down+1/sqrt(2)*self.up)
        else :
            psi0 = np.array([0,1,0])
        resultR,control=self.compute(w,psi0,c,dt)
        ry=[resultR[0][-1],resultR[1][-1],resultR[2][-1]]

        if Periodic_system.qutip:
            psi0=self.up
        else :
            psi0 = np.array([0,0,1])
        resultR,control=self.compute(w,psi0,c,dt)
        rz=[resultR[0][-1],resultR[1][-1],resultR[2][-1]]


        rx,ry,rz=np.array(rx),np.array(ry),np.array(rz)
        M=[rx,ry,rz]
        M=np.array(M).transpose()
        print('M',M)
        return self.get_axis_angle_aux(M)
        
    def get_axis_angle_aux(self,M):
        """compute axis and angle from SO3 matrix M"""
        eigen_values,eigen_vectors=np.linalg.eig(M)
        print('eigen_value',eigen_values)
        h=np.argmin(np.abs(eigen_values-np.ones(3))) # we take the eigen vector associated with the eigen_value 1
        v=eigen_vectors[:,h] # the fixe point
        v=np.real(v) #the imaginary part should be only computational error
        #print(v)
        # now we need the angle
        e1=np.array([1,0,0])
        e2=np.array([0,1,0])
        if np.dot(e1,v)>0.9:
            vector_ortho=np.cross(e2,v)# get an orthogonal vector to v
        else :
            vector_ortho=np.cross(e1,v)
        vector_ortho=vector_ortho/np.linalg.norm(vector_ortho)# we normalize it
        costheta=max(min((np.trace(M)-1)/2,1),-1)
        theta = np.arccos(costheta)
        if np.dot(np.cross(vector_ortho,np.dot(M,vector_ortho)), v)<0 :# the rotation must be between 0 and pi
            v=-v
        #print(theta)
        #print(np.dot(vector_ortho,v))
        #print(np.cross(vector_ortho,np.dot(M,vector_ortho)))
        #print(np.linalg.norm(np.cross(vector_ortho,np.dot(M,vector_ortho))))
        #print(np.cos(theta))
        #print(costheta)
        #testing the result :
        if self.plot_inter==True:
            psi0=sqrt((v[2]+1)/2)*self.up+sqrt(1-(v[2]+1)/2)*np.exp(1j*np.arctan2(v[1],v[0]))*self.down # create the right vector
            resultR=self.compute(w,psi0,c,dt)
            rw=[resultR[0][-1],resultR[1][-1],resultR[2][-1]]
            times=np.arange(0,2*np.pi/w,dt)
            fig = plt.figure()#figsize=(15,22))
            ax = fig.add_subplot(1,2,1)
            plt.title('Evolution real case')
            ax.plot(times, resultR[0], 'r')
            ax.plot(times, resultR[1], 'g')
            ax.plot(times, resultR[2], 'b')
            ax.plot(times, control, 'orange')


            g=np.gradient(resultR[0],dt)/resultR[1]


            #ax.plot(times, g, 'y')
            ax.legend(("sx", "sy", "sz"))

            sphere=Bloch(axes=fig.add_subplot(1,2,2, projection='3d'))
            sphere.add_points([resultR[0],resultR[1],resultR[2]], meth='l')
            sphere.vector_color = ['r']
            #sphere.add_vectors([sin(theta),0,cos(theta)])
            sphere.make_sphere()
            sphere.show()
            plt.show()
        return(v,theta)

    def launch_simu(self,nb_point=40,E=2,dt=0.001,c=1,alpha=0,pts=None):
        if pts is None:
            pts=np.linspace(0,2*pi,nb_point)
        self.E=E
        self.alpha=alpha
        lst_v=[]
        lst_theta=[]
        for tt in pts:
            w=np.abs(2*(E)+self.control_u2(tt))
            print('tt,w',tt,w)
            v,theta=self.get_axis_angle(dt,w,c*self.control_u1(tt))
            #for a better visualisation, we are going to force x>0
            if v[0]>0:
                lst_v.append(v)
                lst_theta.append(theta)
            else :
                lst_v.append(-v)
                lst_theta.append(2*np.pi-theta)#we keep theta in 0:2*pi because we expecte values near pi

        lst_theta=np.array(lst_theta)
        lst_v=np.array(lst_v)
        x=lst_v[:,0]
        y=lst_v[:,1]
        z=lst_v[:,2]
        if self.plot_res:
            self.plot_result(x,y,z,lst_theta)
        return(x,y,z,lst_theta)

    def plot_result(self,x,y,z,lst_theta):
        fig = plt.figure()
        ax = fig.add_subplot(111)#, projection='3d')
        pts=np.linspace(0,2*pi,len(x))
        self.regularize(lst_theta)
        ax.plot(pts, x, 'r')
        ax.plot(pts, y, 'g')
        ax.plot(pts, z, 'b')
        ax.plot(pts, lst_theta, 'orange')
        ax.legend(("sx", "sy", "sz",'angle'))
        plt.show()

    def plot_around_zero(self,ax,x,z,theta):
        # 2D plot (in x,z) giving the rotation angle
        n=len(x)
        lx=np.zeros(n)
        lz=np.zeros(n)
        for k in range(n):
            #print(x[k],z[k],theta[k])
            if theta[k]< np.pi:
                lx[k]=x[k]*theta[k]
                lz[k]=z[k]*theta[k]
            else :
                lx[k]=x[k]*(2*np.pi-theta[k])
                lz[k]=z[k]*(2*np.pi-theta[k])
        ax.plot(lx,lz,'x')
        ax.legend(["axis coordinate (x,z)"])

    def regularize(self,lst_theta):
        """change lst_theta for avoiding jump between 2 pi and 0"""
        i=0
        for k in range(len(lst_theta)-1):
            if np.abs(lst_theta[k+1]-lst_theta[k]+i*2*np.pi)>np.abs(lst_theta[k+1]-lst_theta[k]+2*np.pi+i*2*np.pi):
                i=i+1
            elif np.abs(lst_theta[k+1]-lst_theta[k]+i*2*np.pi)>np.abs(lst_theta[k+1]-lst_theta[k]-2*np.pi+i*2*np.pi):
                i=i-1
            lst_theta[k+1]=lst_theta[k+1]+2*i*np.pi

def monothread():
    test=Periodic_system(False)
    lst_alpha=np.arange(-1,1,0.1)
    res_tot=[]
    #i=int((sys.argv[1]))
    i=10
    alpha=lst_alpha[i]
    (x,y,z,lst_theta)=test.launch_simu(nb_point=1000,E=2,dt=0.0001,c=1,alpha=alpha)
    res_tot.append(x)
    res_tot.append(y)
    res_tot.append(z)
    res_tot.append(lst_theta)
    time.sleep(5*np.random.random(1))#eviter que tout le monde sauvegarde en meme temps
    try:
        path='simu_xyztheta_1000pt-4dt_{:f}'.format(alpha)
        with open(path, 'wb') as fp:
            pickle.dump(res_tot, fp)
        os.system('scp '+path+ ' remi@129.104.219.230:test/')
        os.system('rm '+path) 
        print("done")
    except:
        time.sleep(1)
        os.system('echo "{:d}" >> res/failed3.txt'.format(i))

if __name__=='__main__':
    #monothread()
    test=Periodic_system(False)
    alpha=0
    #pts=np.linspace(0.94,0.99,50)
    #(x,y,z,lst_theta)=test.launch_simu(nb_point=1000,E=2,dt=0.0005,c=1,alpha=alpha,pts=pts)
#    with open('res/simu_xyztheta_zeros_8000pt-4dt_-0.000000','rb') as fg:
    with open('res/simu_xyztheta_zeros_2500_euc_pt-6dt_-0.000000','rb') as fg :
        (x,y,z,lst_theta)=pickle.load(fg)
        print(z)
        x=np.array(x)
        x=x.flatten()
        y=np.array(y)
        y=y.flatten()
        z=np.array(z)
        z=z.flatten()
        lst_theta=np.array(lst_theta)
        lst_theta=lst_theta.flatten()
    fig = plt.figure()
    ax = fig.add_subplot(111)#, projection='3d')
    test.plot_around_zero(ax,x,z,lst_theta)
    plt.show()