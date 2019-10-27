# Python code to simulate buildings in earthquakes
import numpy as np
import matplotlib.pyplot as plt

# Time step
T = 100
dt = 1.5
N = int(T/dt)
t = np.linspace(0,N*dt,N)

# Mass matrix
mf = 1
md = 0.1
M = [[mf, 0, 0], [0, mf, 0], [0, 0, md]]

# Stiffness matrix
kf = 1
kd = 0.1
K = [[2*kf, -kf, 0], [-kf, kf + kd, -kd], [0, -kd, kd]]

# Damping matrix
hf = 0.01
hd = 0.1
H = [[2*hf, -hf, 0], [-hf, hf + hd, -hd], [0, -hd, hd]]

#Force
Force = 1.0

# Initial conditions
x = np.zeros([N,3])
v = np.zeros([N,3])
a = np.zeros([N,3])
x[0],x[1] =[0,0,0],[0,0,0]
v[0],v[1] =[0,0,0],[0,0,0]



def max_amplitude(freq):
    for p in range(2,N):
        F = [Force*np.cos(freq*p*dt),0,0]
        a[p] = np.matmul(np.linalg.inv(M), F - np.matmul(K, x[p-1]) - np.matmul(H, v[p-1])) 
        v[p] = v[p-1]+a[p]*dt
        x[p] = x[p-1]+v[p]*dt
        
    max_amp= [ max(x[int(len(x)*0.8):,floor]) for floor in range(3) ] 
    if abs(freq-1.0)<0.1:
        plt.plot(range(N),x[:,0])
        plt.show()
    return max_amp

freq_sweep = np.linspace(0.1,3,30)
print(freq_sweep)
response = np.array([ max_amplitude(freq) for freq in freq_sweep])

plt.plot(freq_sweep,response[:,0],label = 'floor1')
plt.plot(freq_sweep,response[:,1],label = 'floor2')
plt.plot(freq_sweep,response[:,2],label = 'damper')
plt.legend()
plt.show()











