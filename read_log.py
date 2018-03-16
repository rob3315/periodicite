import numpy as np
import matplotlib.pyplot as plt
x,y=[],[]
fp=open('res/logm.txt')
n=8000
for i in range(n):
    a,b=fp.readline(),fp.readline()
    mini=2
    second=2
    for k in b.split(' '):
    	if len(k)>9:
    		if '\n' in k:
    			k=k.split(']')[0]
    		print(k,i)
    		z=complex(k)
    		a=abs(z-1)
    		if a<mini:
    			second=mini
    			mini=a
    		elif a>second:
    			second=a
    x.append(mini)
    y.append(second)

fig=plt.figure()
plt.plot(np.arange(n),np.log(x))
plt.plot(np.arange(n),np.log(y))
plt.legend(('nearest eigen value to 1 (log |lambda-1|)',"second nearest eigen value to 1 (log |lambda-1|)") )
plt.show()