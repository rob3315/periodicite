import sys
import os
import time
i=int((sys.argv[1]))
with open('todo.txt','r') as f:
    ll=f.readlines()
v=ll[i]
print("run "+v)
time.sleep(i*10)
os.system("python simu_continue.py "+v)
