from qutip import *
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
pi=np.pi
from math import sqrt
import pickle
import os
import time
import sys
class Periodic_system():
    def __init__(self,plot):
        self.up = basis(2, 0)
        self.down=basis(2, 1)
        self.control_u1 = lambda t : 1.5*(1-np.cos(t))
        self.control_u2 = lambda t : -3 *np.cos(t/2)
        self.control_U2= lambda t: -6 * np.sin(t/2) # primitive of u2 such as U2(0)=0
        self.plot_res=plot # plot of the result
        self.plot_inter=False # plot of every configuration
    def compute(self,w,psi0,c,dt):
        H0=(self.E+self.alpha)*sigmaz()
        times=np.arange(0,2*np.pi/w,dt)


        def H1_coeffR(t,args):
            return c*1/2*np.real( np.exp( -1j*(w*t) ))

        #def H1_coeffCx(t,args):
        #    return u1(epsilon*tt) *1/4*np.real( np.exp( -1j*w*t ))

        #def H1_coeffCy(t,args):
        #    return -u1(epsilon*tt) *1/4*np.imag( np.exp( -1j*w*t ))


        HR=[H0,[sigmax(),H1_coeffR]]
        #HC=[H0,[sigmax(),H1_coeffCx],[sigmay(),H1_coeffCy]]

        resultR = mesolve(HR, psi0, times, [], [sigmax(), sigmay(),sigmaz()])
        #= mesolve(HC, psi0, times, [], [sigmax(), sigmay(),sigmaz()])
        #print(H1_coeffR(times,()))
        return(resultR,H1_coeffR(times,()))


    def get_axis_angle(self,dt,w,c):
    # we compute the evolution of ex, ey and ez
    # in fact it is enough to get the evolution of ex and ey because
    #the evolution matrix in in SO_3 

        psi0=(1/sqrt(2)*self.down+1/sqrt(2)*self.up)
        resultR,control=self.compute(w,psi0,c,dt)
        rx=[resultR.expect[0][-1],resultR.expect[1][-1],resultR.expect[2][-1]]

        psi0=(1j/sqrt(2)*self.down+1/sqrt(2)*self.up)
        resultR,control=self.compute(w,psi0,c,dt)
        ry=[resultR.expect[0][-1],resultR.expect[1][-1],resultR.expect[2][-1]]


        psi0=self.up
        resultR,control=self.compute(w,psi0,c,dt)
        rz=[resultR.expect[0][-1],resultR.expect[1][-1],resultR.expect[2][-1]]


        rx,ry,rz=np.array(rx),np.array(ry),np.array(rz)
        M=[rx,ry,rz]
        M=np.array(M).transpose()
        eigen_values,eigen_vectors=np.linalg.eig(M)
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
        costheta=max(min((1-np.trace(M))/2,1),-1)
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
            resultR,control=self.compute(w,psi0,c,dt)
            rw=[resultR.expect[0][-1],resultR.expect[1][-1],resultR.expect[2][-1]]
            times=np.arange(0,2*np.pi/w,dt)
            fig = plt.figure()#figsize=(15,22))
            ax = fig.add_subplot(1,2,1)
            plt.title('Evolution real case')
            ax.plot(times, resultR.expect[0], 'r')
            ax.plot(times, resultR.expect[1], 'g')
            ax.plot(times, resultR.expect[2], 'b')
            ax.plot(times, control, 'orange')


            g=np.gradient(resultR.expect[0],dt)/resultR.expect[1]


            #ax.plot(times, g, 'y')
            ax.legend(("sx", "sy", "sz"))

            sphere=Bloch(axes=fig.add_subplot(1,2,2, projection='3d'))
            sphere.add_points([resultR.expect[0],resultR.expect[1],resultR.expect[2]], meth='l')
            sphere.vector_color = ['r']
            #sphere.add_vectors([sin(theta),0,cos(theta)])
            sphere.make_sphere()
            sphere.show()
            plt.show()
        return(v,theta)

    def launch_simu(self,nb_point=40,E=2,dt=0.001,c=1,alpha=0):
    
        pts=np.linspace(0,2*pi,nb_point)
        self.E=E
        self.alpha=alpha
        lst_v=[]
        lst_theta=[]
        for tt in pts:
            w=np.abs(2*(E)+self.control_u2(tt))
            print(tt,w)
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

    def regularize(self,lst_theta):
        """change lst_theta for avoiding jump between 2 pi and 0"""
        i=0
        for k in range(len(lst_theta)-1):
            if np.abs(lst_theta[k+1]-lst_theta[k]+i*2*np.pi)>np.abs(lst_theta[k+1]-lst_theta[k]+2*np.pi+i*2*np.pi):
                i=i+1
            elif np.abs(lst_theta[k+1]-lst_theta[k]+i*2*np.pi)>np.abs(lst_theta[k+1]-lst_theta[k]-2*np.pi+i*2*np.pi):
                i=i-1
            lst_theta[k+1]=lst_theta[k+1]+2*i*np.pi

if __name__=='__main__':
    test=Periodic_system(False)
    lst_alpha=np.arange(-1,1,0.1)
    res_tot=[]
    i=int((sys.argv[1]))
    alpha=lst_alpha[i]
    (x,y,z,lst_theta)=test.launch_simu(nb_point=1000,E=2,dt=0.0001,c=1,alpha=alpha)
    res_tot.append(x)
    res_tot.append(y)
    res_tot.append(z)
    res_tot.append(lst_theta)
    time.sleep(5*np.random.random(1))#eviter que tout le monde sauvegarde en meme temps
    try:
        with open('simu_xyztheta_1000pt_-4dt_{:f}'.format(alpha), 'wb') as fp:
            pickle.dump(res_tot, fp)
        os.system('scp '+path+ ' remi@129.104.219.230:test/')
        os.system('rm '+path) 
        print("done")
    except:
        time.sleep(1)
        os.system('echo "{:d}" >> res/failed2.txt'.format(i))