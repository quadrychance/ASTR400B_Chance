# Homework 4
# Center of Mass Position and Velocity

# import modules
import numpy as np
import astropy.units as u
from ReadFile import Read
from astropy.table import Table
from astropy.io import ascii
class CenterOfMass:

    def __init__(self, filename, ptype):
        # read in the file and particle type
        self.time, self.total, self.data = Read(filename)

        #create an array to store indexes of particles of desired Ptype
        self.index = np.where(self.data['type'] == ptype)

        # store the mass, positions, velocities of only the particles of the given type
        self.m = self.data['m'][self.index]

        ##### PLACE other particle properties here: x,y,z,vx,vy,vz #####


                                                                    #initializing other things that go into the class
        self.x = self.data['x'][self.index]#*u.kpc                  # units broke something with how variable were getting read so i just took them out
        self.y = self.data['y'][self.index]#*u.kpc
        self.z = self.data['z'][self.index]#*u.kpc

        self.vx = self.data['vx'][self.index]#*u.km/u.s
        self.vy = self.data['vy'][self.index]#*u.km/u.s
        self.vz = self.data['vz'][self.index]#*u.km/u.s



        #return  xnew, ynew, znew, vxnew, vynew, vznew,
    def total_mass(self):
        #Note: you can add other keyword arguments into the function, but 'self' must be first
        return np.sum(self.m)*u.Msun*1e10

    ##### PLACE OTHER FUNCTIONS BELOW #####

    def COMdefine(self, i,j,k,m):                             #general equation for COM
        xcom = np.sum(i*m)/(np.sum(self.m))
        ycom = np.sum(j*m)/(np.sum(self.m))
        zcom = np.sum(k*m)/(np.sum(self.m))
        return xcom,ycom,zcom
        print(xcom)
    def COM_P(self,delta):

        #com_pos = self.COMdefine(self.x,self.y,self.z,self.m)

        xcom, ycom, zcom = self.COMdefine(self.x,self.y,self.z,self.m) #putting the particles into the com frame
        rcom = np.sqrt(xcom**2+ycom**2+zcom**2)
        xnew = self.x-xcom
        ynew = self.y-ycom
        znew = self.z-zcom
        rnew = np.sqrt((xnew)**2+(ynew)**2+(znew)**2)  #com r magnitude

        difference = 1000                                  #the function needs a starting place, so i chose 1000 kpc
        delta = 10                                         #delta is the codition for convergence

        rmax = np.max(rnew)/2
        while (difference>delta):                          ##the while loop is to iterate through distances down from the starting distance





            index = np.where(rnew < rmax )                #this index will cut out all particles outside the last radius the function used

            xnew = self.x[index]
            ynew = self.y[index]
            znew = self.z[index]
            mnew = self.m[index]

            xcom2, ycom2, zcom2 = self.COMdefine(xnew, ynew, znew, mnew) #xcom2 etc... are where the next direction vector is stored until i replace xcom with it.
            rcom2 = np.sqrt((xcom2**2)+(ycom2**2)+(zcom2**2))
            difference = np.absolute(rcom2-rcom)
            xnew2 = self.x-xcom2
            ynew2 = self.y-ycom2
            znew2 = self.z-zcom2
            rnew = np.sqrt((xnew2**2)+(ynew2**2)+(znew2**2))

            rmax = rmax/2
            rcom = rcom2
            xcom = xcom2
            ycom = ycom2
            zcom = zcom2


            return(xcom, ycom, zcom)                                   #this continues until the change from xcom2 to xcom is delta

            print(delta)

    def COM_V(self):
        xcom, ycom, zcom = self.COMdefine(self.x,self.y,self.z,self.m)  #now that i have the com position i can get the velocity
        rcom = np.sqrt(xcom**2+ycom**2+zcom**2)
        xnew = self.x-xcom
        ynew = self.y-ycom
        znew = self.z-zcom
        rnew = np.sqrt((xnew)**2+(ynew)**2+(znew)**2)


        index = np.where(rnew <= 15 )

        vxnew = self.vx[index]
        vynew = self.vy[index]
        vznew = self.vz[index]
        mnew = self.m[index]
        vxcom,vycom,vzcom = self.COMdefine(self.vx, self.vy, self.vz, self.m)
        return(vxcom, vycom, vzcom)



# EXAMPLE OF USING A CLASS
##########################

# Create a Center of mass object for the galaxies
MWCOM = CenterOfMass("MW_000.txt", 2)
M31COM = CenterOfMass("M31_000.txt", 2)
M33COM = CenterOfMass("M33_000.txt", 2)

# Calculate quantities for MW data
MW_mass = MWCOM.total_mass()
MW_pos = MWCOM.COM_P(10)
MW_vel = MWCOM.COM_V()

M31_mass = M31COM.total_mass()  #M31 data
M31_pos = M31COM.COM_P(10)
M31_vel = M31COM.COM_V()

M33_mass = M33COM.total_mass() #M33 data
M33_pos = M33COM.COM_P(10)
M33_vel = M33COM.COM_V()

pos = [MW_pos, M31_pos, M33_pos]   #astropy.tables likes data in lists
mass = [MW_mass, M31_mass, M33_mass]
vel = [MW_vel, M31_vel, M33_vel]

t = Table([pos, mass, vel], names=('Galaxy Position', 'Galaxy Mass', 'Galaxy Velocity'), dtype = [('object'), ('object'),('object')],meta={'name': 'Mass in trillion mSun'}) #write class outputs to table
print(t)

ascii.write(t, 'table.dat', format = 'latex', overwrite = True)

print(np.sqrt(np.asarray(MWCOM.COM_P(10))**2)-np.asarray(M31COM.COM_P(10))**2)  # these return the squared components that i just added and took the square root of by hand because im tired

print(np.sqrt(np.asarray(MWCOM.COM_V())**2)-np.asarray(M31COM.COM_V())**2)

print(np.sqrt(np.asarray(M31COM.COM_P(10))**2)-np.asarray(M33COM.COM_P(10))**2)
print(np.sqrt(np.asarray(M31COM.COM_V())**2)-np.asarray(M33COM.COM_V())**2)
