# Python code to simulate buildings in earthquakes
import numpy as np
import matplotlib.pyplot as plt

# Mode
mode = "imp" # imp, freq or damp
integrator = "euler" # euler or verlet

# Time step
T = 40
dt = 0.01
N = int(T/dt)
t = np.linspace(0, N*dt, N)

# Mass matrix
mf = 3.94
md = 0.1
M = [[mf, 0, 0], [0, mf, 0], [0, 0, md]]

# Stiffness matrix
kf = 2100
kd = 865
K = [[2*kf, -kf, 0], [-kf, kf + kd, -kd], [0, -kd, kd]]

# Damping matrix
hf = 2.5
hd = 0.5

#Force
Force = 1.0

#Undamped natural frequency
resonance = 14.8

# Initial conditions
x = np.zeros([N,3])
v = np.zeros([N,3])
a = np.zeros([N,3])
x[0],v[0] =[0, 0, 0], [0, 0, 0]

def run(freq=0, damp=hd):
    
    H = [[2*hf, -hf, 0], [-hf, hf + damp, -damp], [0, -damp, damp]]
    
    for p in range(1, N):
        if mode == "freq" or mode == "damp":
            F = [Force * np.cos(freq * p * dt), 0, 0]
        elif mode == "imp" and p == 1:
            F = [1/dt, 0, 0]
        else:
            F = [0, 0, 0]
    
        if integrator == "verlet":
            a[p] = np.matmul(np.linalg.inv(M), F - np.matmul(K, x[p-1]) - np.matmul(H, v[p-1])) 
            if p != 1:
                x[p] = 2 * x[p-1] - x[p-2] + dt * dt * a[p]
            else:
                x[p] = 2 * x[p-1] + dt * dt * a[p]
            v[p] = (x[p] - x[p-1])/dt

        elif integrator == "euler":
            a[p] = np.matmul(np.linalg.inv(M), F - np.matmul(K, x[p-1]) - np.matmul(H, v[p-1])) 
            v[p] = v[p-1] + a[p]*dt
            x[p] = x[p-1] + v[p]*dt


def max_amplitude(freq, damp):
    run(freq, damp)
    max_amp = [ max(x[int(len(x)*0.8):,floor]) for floor in range(3) ]
    
    return max_amp

# Formating

font = {'weight' : 'normal',
        'size'   : 22}

plt.rc('font', **font)

if mode == "freq":
    freq_sweep = np.linspace(0.1, 60, 50)
    response = np.array([max_amplitude(freq, hd) for freq in freq_sweep])
    lines = plt.plot(freq_sweep,response)

    plt.legend(lines, ('floor1', 'floor2', 'damper'))
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Amplitude (m)')

    plt.show()
    
elif mode == "imp":
    run()
    lines = plt.plot(t, x)

    plt.legend(lines, ('floor1', 'floor2', 'damper'))
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude (m)')

    plt.show()

if mode == "damp":
    damp_sweep = np.linspace(0.1, 15, 10)
    response = np.array([max_amplitude(resonance, damp) for damp in damp_sweep])
    lines = plt.plot(damp_sweep,response)

    plt.legend(lines, ('floor1', 'floor2', 'damper'))
    plt.xlabel('Lambda (Ns/m)')
    plt.ylabel('Amplitude (m)')

    plt.show()




