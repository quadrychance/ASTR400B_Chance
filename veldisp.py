import numpy as np
import astropy.units as u
from astropy.constants import G
from ReadFile import Read
from CenterOfMass import CenterOfMass
from MassDist import MassProfile
 #from lab7
 # Create a COM of object for M31 Disk Using Code from Assignment 4
# modified CenterOfMass so that it doesn't return rounded values
COMD = CenterOfMass("M31_000.txt",2)

# Compute COM of M31 using disk particles
COMP = COMD.COM_P(0.1,4.0)
COMV = COMD.COM_V(COMP[0],COMP[1],COMP[2])

# Determine positions of disk particles relative to COM
xD = COMD.x - float(COMP[0]/u.kpc)
yD = COMD.y - float(COMP[1]/u.kpc)
zD = COMD.z - float(COMP[2]/u.kpc)

# total magnitude
rtot = np.sqrt(xD**2 + yD**2 + zD**2)

# Determine velocities of disk particles relatiev to COM motion
vxD = COMD.vx - float(COMV[0]/u.km*u.s)
vyD = COMD.vy - float(COMV[1]/u.km*u.s)
vzD = COMD.vz - float(COMV[2]/u.km*u.s)

# total velocity
vtot = np.sqrt(vxD**2 + vyD**2 + vzD**2)

# Vectors
r = np.array([xD,yD,zD]).T
v = np.array([vxD,vyD,vzD]).T

index = np.where(np.abs(vtot) < 300)

def RotateFrame(posI,velI):
    # input:  3D array of positions and velocities
    # returns: 3D array of rotated positions and velocities such that j is in z direction

    # compute the angular momentum
    L = np.sum(np.cross(posI,velI), axis=0)
    # normalize the vector
    L_norm = L/np.sqrt(np.sum(L**2))


    # Set up rotation matrix to map L_norm to z unit vector (disk in xy-plane)

    # z unit vector
    z_norm = np.array([0, 0, 1])

    # cross product between L and z
    vv = np.cross(L_norm, z_norm)
    s = np.sqrt(np.sum(vv**2))

    # dot product between L and z
    c = np.dot(L_norm, z_norm)

    # rotation matrix
    I = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    v_x = np.array([[0, -vv[2], vv[1]], [vv[2], 0, -vv[0]], [-vv[1], vv[0], 0]])
    R = I + v_x + np.dot(v_x, v_x)*(1 - c)/s**2

    # Rotate coordinate system
    pos = np.dot(R, posI.T).T
    vel = np.dot(R, velI.T).T

    return pos, vel
r_new, v_new =  RotateFrame(r[index],v[index])

sigma_x = np.mean(np.sqrt(np.abs(v_new[:,0]**2)))
sigma_y = np.mean( np.sqrt(np.abs(v_new[:,1]**2)))
sigma_z = np.mean(np.sqrt(np.abs(v_new[:,2]**2)))

#get circular of disk particles velcoity

MW_mp = np.zeros(10, 0, 800)
M31_mp = np.zeros(10, 0, 800)

MW_CV = np.zeros(1,0,50)
M31_CV = np.zeros(1,0,50)

for n in range(0,800):           #trying to get mp for each snapshot
    MW_mp = MassProfile('MW', n)
    M31_mp = MassProfile('M31', n)
    MW_CV = MW_mp.CircularVelocity(2,10 )
    M31_CV = M31_mp.CircularVelocity(2, 10)
print(MW_CV)
