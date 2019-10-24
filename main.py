# Python code to simulate buildings in earthquakes
import numpy as np

# Time step

dt = 0.1

# Mass matrix
m1 = 1
m2 = 1
md = 1

M = np.matrix([[m1, 0, 0], [0, m2, 0], [0, 0, md]])

# Stiffness matrix

k1 = 1
k2 = 1
kd = 1

K = np.matrix([[k1 + k2, -k2, 0], [-k1, k1 + kd, -kd], [0, -kd, kd]])

# Damping matrix

h1 = 1
h2 = 1
hd = 1

H = np.matrix([[h1 + h2, -h2, 0], [-h1, h1 + hd, -hd], [0, -hd, hd]])

# Initial conditions

x1 = 0
x2 = 0
xd = 0

v1 = 1
v2 = 0
vd = 0

x0 = np.matrix([x1],[x2],[xd])
v0 = np.matrix([v1],[v2],[vd])

x = [x0]
v = [v0]

for p in range(2,N):
    a[p] = np.dot(np.linalg.inv(M), np.dot(K, x[p-1]) + np.dot(H, v[p-1])) 
    x[p] = 2 * x[p-1] - x[p-2] + dt * dt * a[p]
    v[p] = (x[p] - x[p-1])/dt

print(x)
