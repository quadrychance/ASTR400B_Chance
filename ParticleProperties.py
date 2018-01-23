import numpy as np
import astropy.units as u
from ReadFile import Read

#def ParticleProperties(pid): <- making it a function broke my program so i changed it back
t,p,d = (Read('MW_000.txt'))

pid = 99 # particle number


#get 3d distance and velocity. i'm unsure how to get direction
#using 3d distance as a variable doesnt work?
A = ((d['x'][pid]**2)+(d['y'][pid]**2)+(d['z'][pid]**2))**(.5)*1000*u.pc
B = ((d['vx'][pid]**2)+(d['vy'][pid]**2)+(d['vz'][pid]**2))**(.5)

#Mass in m_solar
C = ((d['m'][pid])*(10**10))

print('3D distance = ')
#convert kpc to lyr
print(np.around(A.to(u.lyr),3))

print('3D velocity = ')
print(np.around(B,3))
print('mass in solar masses')
print(C)
