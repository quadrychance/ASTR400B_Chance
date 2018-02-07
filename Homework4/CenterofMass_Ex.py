# Homework 4 
# Center of Mass Position and Velocity

# import modules
import numpy as np
import astropy.units as u
from ReadFile import Read


class CenterOfMass:
   
    def __init__(self, filename, ptype):
        # read in the file and particle type
        self.time, self.total, self.data = Read(filename)
            
        #create an array to store indexes of particles of desired Ptype
        self.index = np.where(self.data['type'] == ptype)
    
        # store the mass, positions, velocities of only the particles of the given type
        self.m = self.data['m'][self.index]
    
        ##### PLACE other particle properties here: x,y,z,vx,vy,vz #####
      
    def total_mass(self):
        #Note: you can add other keyword arguments into the function, but 'self' must be first
        return np.sum(self.m)*u.Msun*1e10
        
    ##### PLACE OTHER FUNCTIONS BELOW #####





# EXAMPLE OF USING A CLASS
##########################

# Create a Center of mass object for the MW
MWCOM = CenterOfMass("MW_000.txt", 2)

# Calculate quantities for MW data
MW_mass = MWCOM.total_mass()
print("MW Disk Mass:", MW_mass)
