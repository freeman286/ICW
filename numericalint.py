import numpy as np
import matplotlib.pyplot as plt

# Python code to simulate buildings in earthquakes
import numpy as np

# Time step
T = 50
dt = 0.001
N = int(T/dt)
# Mass matrix
m1 = 1
m2 = 1
md = 1

M = [[m1, 0, 0], [0, m2, 0], [0, 0, md]]

# Stiffness matrix

k1 = 1
k2 = 1
kd = 1

K = [[k1 + k2, -k2, 0], [-k1, k1 + kd, -kd], [0, -kd, kd]]

# Damping matrix

h1 = 1
h2 = 1
hd = 1

H = [[h1 + h2, -h2, 0], [-h1, h1 + hd, -hd], [0, -hd, hd]]

# Initial conditions

x1 = 0
x2 = 0
xd = 0

v1 = 1
v2 = 0
vd = 0

x0 =[x1,x2,xd]
v0 =[v1,v2,vd]

x = np.zeros([N,3])
v = np.zeros([N,3])
a = np.zeros([N,3])


x[0] = x0
v[0] = v0
x[1] = x[0]+v[0]*dt
v[1] = v0

print(np.matmul(H,v0))


for p in range(2,N):
    a[p] = -np.matmul(np.linalg.inv(M), np.matmul(K, x[p-1]) + np.matmul(H, v[p-1])) 
    x[p] = 2 * x[p-1] - x[p-2] + dt * dt * a[p]
    v[p] = (x[p] - x[p-1])/dt

print(x[:,0])

print(len(x[:,0]),x[:,0])
plt.plot(range(len(x[:,0])),x[:,0])
plt.show()

T = 50
dt = 0.01
N = int(T/dt)
t = np.linspace(0,50,num = N)
x = [np.cos(t),np.sin(2*t),np.sin(5*t)]
print("Max Frequency x0:")
FT = np.abs(np.fft.fft(x[0]))
print(FT)
#plt.plot(range(len(FT)),FT)
#plt.show()
print(max(np.fft.fft(x[0])))