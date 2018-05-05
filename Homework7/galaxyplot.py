import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

MW = np.genfromtxt('Orbits_MW.txt', names = True )
M31 = np.genfromtxt('Orbits_M31.txt', names = True)
M33 = np.genfromtxt('Orbits_M33.txt', names = True)

MW_M31V = np.sqrt((MW['vx']-M31['vx'])**2 + (MW['vy']-M31['vy'])**2 + (MW['vz']-M31['vz'])**2)
MW_M31R = np.sqrt((MW['x']-M31['x'])**2 + (MW['y']-M31['y'])**2 + (MW['z']-M31['z'])**2)

plt.plot(MW['t'], MW_M31R)
plt.xlabel('Time')
plt.ylabel('Distance')
plt.title('MW & M31')

fig = plt.figure(figsize=(10,10))
#ax = plt.subplot(121)
plt.plot(MW['t'], MW_M31V)
plt.xlabel('Time')
plt.ylabel('Velocity')
plt.title('MW & M31')




M31_M33V = np.sqrt((M31['vx']-M33['vx'])**2 + (M31['vy']-M33['vy'])**2 + (M31['vz']-M33['vz'])**2)
M31_M33R =  np.sqrt((M31['x']-M33['x'])**2 + (M31['y']-M33['y'])**2 + (M31['z']-M33['z'])**2)


fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

plt.plot(M31['t'], M31_M33R)
plt.xlabel('Time')
plt.ylabel('Distance')
plt.title('M33 & M31')

fig = plt.figure(figsize=(10,10))
#ax = plt.subplot(121)

plt.plot(M31['t'], M31_M33V)
plt.xlabel('Time')
plt.ylabel('Velocity')
plt.title('M33 & M31')
plt.show()
