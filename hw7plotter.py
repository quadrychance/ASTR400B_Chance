import matplotlib.pyplot as plt
from M33_Analytic import M33AnalyticOrbit
import numpy as np


M31Orbit = np.genfromtxt('Orbits_M31.txt',names=True)
M33Orbit = np.genfromtxt('Orbits_M33.txt',names=True)
M33Analytic = np.genfromtxt('M33Analytic',names=True)


sim_sep = np.sqrt((M33Orbit['x']-M31Orbit['x'])**2+(M33Orbit['y']-M31Orbit['y'])**2+(M33Orbit['z']-M31Orbit['z'])**2)

a_sep = np.sqrt(M33Analytic['x']**2+M33Analytic['y']**2+M33Analytic['z']**2)

plt.xkcd()
plt.figure(1,figsize=(10,10))
plt.plot(M31Orbit['t']/10, sim_sep, c= 'r', label ='Sim')  #had some weird timescale issue here
plt.plot(M33Analytic['t'], a_sep, c='b', label ='Analytic')
plt.title('Seperation of M33 and M31')
plt.xlim(0,20)

plt.legend()
plt.show()

sim_v = np.sqrt((M31Orbit['vx']-M33Orbit['vx'])**2+(M31Orbit['vy']-M33Orbit['vy'])**2+(M31Orbit['vz']-M33Orbit['vz'])**2)
a_v = np.sqrt(M33Analytic['vx']**2+M33Analytic['vy']**2+M33Analytic['vz']**2)
plt.figure(1,figsize=(10,10))
plt.plot(M31Orbit['t']/10, sim_v, c= 'r', label='Sim')
plt.plot(M33Analytic['t'], a_v, c='b', label='Analytic')
plt.title('Velocity Dispersion between M33 and M31')
plt.legend()

plt.show()
