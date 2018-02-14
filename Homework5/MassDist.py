#THIS DOES NOT RUN.
import numpy as np
import astropy.units as u
from ReadFile import Read
from astropy.table import Table
from astropy.io import ascii
from astropy.constants import G
from CenterofMass import CenterOfMass
from matplotlib import pyplot as plt



radii = np.linspace(.01, 30, 100)
class MassProfile:

    def __init__(self, galaxy, snap):

        ilbl = '000' + str(snap)       #FILENAME RECONSTRUCTION

        ilbl = ilbl[-3:]
        self.filename = '%s_'%(galaxy) + ilbl + '.txt'

         # read in the file and particle type
        self.time, self.total, self.data = Read(self.filename)




        self.gname = galaxy

    def MassEnclosed(self, ptype, radii):   #takes in particle type and the radius array and returns mass enclosed within r

        self.index = np.where(self.data['type'] == ptype)
        self.m = self.data['m'][self.index]
        #print(self.m)

        self.x = self.data['x'][self.index]#*u.kpc   #reading in data
        self.y = self.data['y'][self.index]#*u.kpc
        self.z = self.data['z'][self.index]#*u.kpc



        COM1 = CenterOfMass(self.filename, ptype)

        R = COM1.COM_P(10)#*u.kpc
        xnew = self.x - R[0]
        ynew = self.y - R[1]
        znew = self.z - R[2]

        mass_enclosed = np.zeros(len(radii))
        rnew = np.sqrt(xnew**2+ynew**2+znew**2)

        for i in range(np.size(radii)):
            print(radii[i],rnew)
            r_enc = np.where (rnew > radii[i])#*u.kpc)

            mass_enclosed[i] = np.sum(self.m[r_enc])

            print(np.sum(self.m[r_enc]), i, 'look')

            return(mass_enclosed)#*u.Msun)




    def MassEnclosedTotal(self, radii):             #adds upt MassEnclosed for all the particle types
        return MassEnclosed(1, radii) + MassEnclosed(2, radii) + MassEnclosed(3, radii)




    def HernquistMass(self, radii, a, Mhalo):   #function takes in radius array, scale facto, and halo massr returns the hernquist mass profile


        return (Mhalo*radii**2)/(a + radii)**2

    def CircularVelocity(self, ptype, radii):     #TAKES IN PARTICLE TYPE AND RADIUS ARRAY TO RETURN THE CIRCULART ORBITAL VELOCITY OF PARTICLES AT A RADIUS IN RADII


        return(sqrt(MassEnclosedTotal(radii)*G/radii))


    def HernquistVCirc(self, radii, a, Mhalo): #USES HERQUIST MASS INSTEAD OF THE ENCLOSED MASS TO DO THE SAME


        return((Mhalo*radii**2)/(a + radii)**2*G/radii**2)

r_mp = np.linspace(.1,30,100)

MW_mp = MassProfile('MW', 000)
M31_mp = MassProfile('M31', 000)
M33_mp = MassProfile('M33', 000)

#MASS PROFILE PLOTS
plt.subplot(211)
plt.plot(r_mp, MW_mp.MassEnclosed(1, radii), c = 'k')
plt.plot(r_mp, MW_mp.MassEnclosed(2, radii), c = 'b')
plt.plot(r_mp, MW_mp.MassEnclosed(3, radii), c = 'r')
plt.plot(r_mp, MW_mp.HernquistMass(radii, 1, 1.97E12), c = 'y')
plt.plot(r_mp, MW_mp.MassEnclosedT(radii), c = 'g', marker = '^')
plt.yscale('log')
plt.title('MW')

plt.subplot(212)
plt.plot(r_mp, M31_mp.MassEnclosed(1, radii), c = 'k')
plt.plot(r_mp, M31_mp.MassEnclosed(2, radii), c = 'b')
plt.plot(r_mp, M31_mp.MassEnclosed(3, radii), c = 'r')
plt.plot(r_mp, M31_mp.HernquistMass(radii, 1, 1.97E12), c = 'y')
plt.plot(r_mp, M31_mp.MassEnclosedT(radii), c = 'g', marker = '^')
plt.yscale('log')
plt.title('M31')

plt.subplot(222)
plt.plot(r_mp, M33_mp.MassEnclosed(1, radii), c = 'k')
plt.plot(r_mp, M33_mp.MassEnclosed(2, radii), c = 'b')
plt.plot(r_mp, M33_mp.MassEnclosed(3, radii), c = 'r')
plt.plot(r_mp, M33_mp.HernquistMass(radii, 1, 1.97E12), c = 'y')
plt.plot(r_mp, M33_mp.MassEnclosedT(radii), c = 'g', marker = '^')

plt.yscale('log')
plt.title('M33')
plt.xlabel(radius (kpc))


plt.yscale('log')
plt.show()

#ROTATION CURVE PLOTS
plt.subplot(211)
plt.plot(r_mp, MW_mp.CircularVelocity(1, radii), c = 'k')
plt.plot(r_mp, MW_mp.CircularVelocity(2, radii), c = 'b')
plt.plot(r_mp, MW_mp.CircularVelocity(3, radii), c = 'r')
plt.plot(r_mp, MW_mp.HernquistVCircy(radii, 1, 1.97))

plt.yscale('log')
plt.title('MW')

plt.subplot(212)
plt.plot(r_mp, M31_mp.CircularVelocity(1, radii), c = 'k')
plt.plot(r_mp, M31_mp.CircularVelocity(2, radii), c = 'b')
plt.plot(r_mp, M31_mp.CircularVelocity(3, radii), c = 'r')
plt.plot(r_mp, M31_mp.HernquistVCircy(radii, 1, 1.92))

plt.yscale('log')
plt.title('M31')

plt.subplot(222)
plt.plot(r_mp, M33_mp.CircularVelocity(1, radii), c = 'k')
plt.plot(r_mp, M33_mp.CircularVelocity(2, radii), c = 'b')
plt.plot(r_mp, M33_mp.CircularVelocity(3, radii), c = 'r')
plt.plot(r_mp, M33_mp.HernquistVCircy(radii, 1, 0.186))

plt.yscale('log')
plt.title('M33')
plt.xlabel(radius (kpc))

plt.show()
