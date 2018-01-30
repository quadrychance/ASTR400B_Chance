#Quadry Chance
#Comparing mass of galaxy components


import numpy as np
import astropy.units as u
from ReadFile import Read
from astropy.table import Table
from astropy.io import ascii

def ComponentMass(filename, p_type):
    t,p,d = (Read(filename))   # use readfile function


    index = np.where(d['type'] == p_type) #index the file according to particle type
    m_comp = d['m'][index]*1e-2  #get the numbers to look right
    sum_m_comp = np.sum(m_comp)  #sum all components of the same index
    return(sum_m_comp)
#here i'm running this function on all the in files. Ideally I wouldn't hardcode this, but that more complicated than its worth I think
halo1 = ComponentMass('MW_000.txt',1)
halo2 = ComponentMass('M31_000.txt',1)
halo3 = ComponentMass('M33_000.txt',1)

disk1 = ComponentMass('MW_000.txt',2)
disk2 = ComponentMass('M31_000.txt',2)
disk3 = ComponentMass('M33_000.txt',2)

bulge1 = ComponentMass('MW_000.txt',3)
bulge2 = ComponentMass('M31_000.txt',3)
bulge3 = ComponentMass('M33_000.txt',3)

#place the component masses into a list to be able to use tables
halo_masses = [halo1, halo2, halo3]
disk_masses = [disk1, disk2, disk3]
bulge_masses = [bulge1, bulge2, bulge3]

#add up componet masses for total mass
total_mass1 = halo_masses[0] + disk_masses[0] + bulge_masses[0]
total_mass2 = halo_masses[1] + disk_masses[1] + bulge_masses[1]
total_mass3 = halo_masses[2] + disk_masses[2] + bulge_masses[2]
#stick them into another list
total_masses = [total_mass1, total_mass2, total_mass3]

#galaxy name list
galaxy_names = ['Milky Way', 'M31', 'M33']

#baryonic mass fraction calculation
bar_mass1 = (bulge_masses[0] + disk_masses[0])/total_mass1
bar_mass2 = (bulge_masses[1] + disk_masses[1])/total_mass2
bar_mass3 = (bulge_masses[2] + disk_masses[2])/total_mass3

bar_masses = [bar_mass1, bar_mass2, bar_mass3]

#print(halo_masses)


#astropy.table function. this made making the table super easy. it doesnt like arrays in other arrys, so i had to set the to be read as objects.
t = Table([galaxy_names, halo_masses, disk_masses, bulge_masses, bar_masses, total_masses], names=('Galaxy Name', 'Halo Mass', 'Disk Mass', 'Bulge Mass','Baryon Fraction', 'Total Mass'), dtype = [('object'), ('object'),('object'),('object'),('object'),('object')],meta={'name': 'Mass in trillion mSun'})
print(t)

#another useful astropy thing. you can use the ascii.write function to write to a lot of different types of files. this ouputs a .dat file with the LaTeX code in it.
ascii.write(t, 'table.dat', format = 'latex', overwrite = True)
