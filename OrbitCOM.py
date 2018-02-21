import numpy as np
import astropy.units as u
from ReadFile import Read
from astropy.table import Table
from astropy.io import ascii
from astropy.constants import G
from CenterOfMass import CenterOfMass
from matplotlib import pyplot as plt

def OrbitCOM(galaxy, start, end, n):
    fileout = 'Orbits_M31.txt'
    shape = (int(end/n)+1, 7)   #i think just putting row,column for the shape makes np.zeros think column is the dtype
    Orbit = np.zeros(shape)
    delta = 1
    VolDec = 5


    for i in np.arange(start, end+n, n):
        ilbl = '000' + str(i)
        # remove all but the last 3 digits
        ilbl = ilbl[-3:]
        # create filenames
        filename='%s_'%(galaxy) + ilbl + '.txt'
        time, total, data = Read(filename)

        COM = CenterOfMass(filename,2)

        time = float(COM.time/u.Myr)
        #x = float(COM.x/u.kpc)
        #print(int(i/n))

        R = COM.COM_P(delta,VolDec)
        V = COM.COM_V(R[0],R[1],R[2])
        Orbit[int(i/n),0] = time/1000

        Orbit[int(i/n),1] = float(R[0]/u.kpc)
        Orbit[int(i/n),2] = float(R[1]/u.kpc)
        Orbit[int(i/n),3] = float(R[2]/u.kpc)


        Orbit[int(i/n),4] = float(V[0]/u.km*u.s)
        Orbit[int(i/n),5] = float(V[1]/u.km*u.s)
        Orbit[int(i/n),6] = float(V[2]/u.km*u.s)
        print(i)
        print(Orbit)
    np.savetxt(fileout, Orbit, header='t x y z vx vy vz',comments='#', fmt=['%.2f', '%.2f' ,'%.2f','%.2f','%.2f','%.2f ','%.2f '])


OrbitCOM('M31',0, 345, 5)
