import numpy as np
import astropy.units as u


def Read(filename):
#open file
    file = open(filename, 'r')
#read file
    line1 = file.readline()
#split line 1 on whitespace
    label, value =line1.split()
    time = float(value)*10*u.Myr
#split line 2 on whitespace
    line2 = file.readline()
    label, value = line2.split()
    p_number = int(value)
#close file
    file.close()
#generate np array of datatable, skipping the header and keeping the column names rather than using regular indexing
    data = np.genfromtxt(filename, dtype = None, names = True, skip_header = 3)


    return(time,p_number,data)
