import numpy as np
import astropy.units as u
from ReadFile import Read

class M33AnalyticOrbit:

    def __init__(self,filename):


        self.G = 4.498768e-6 #define constant

        # x, y, z postion of M33 from M31 frame
        self.x = -98
        self.y = -120
        self.z = -127
        #vx, vy, vz of M33 from M31 frame
        self.vx = -28
        self.vy = 174
        self.vz = 92

        #M31 data
        self.rd = 5
        self.Mdisk = 0.12e12
        self.rbulge = 1
        self.Mbulge = 0.019e12
        self.rhalo = 62
        self.Mhalo = 1.921e12

    def HernquistAccel(self,M,r_a,x,y,z,acomp):  #input bulge+halo mass, scale length, positions, wanted accelerationcomponent
        if acomp == 'z':
            c = z
        if acomp == 'y':
            c = y
        if acomp == 'x':
            c = x
        r = np.sqrt(x**2+y**2+z**2)         #position vector
        a = -(self.G*M*c)/(r*(r_a+r)**2)    #hernquist acceleration


        return a

    def MiyamotoNagaiAccel(self,M,rd,x,y,z,acomp):

        H = rd/5           #disk scale height
        B = rd+np.sqrt(z**2+H**2)
        R = np.sqrt(x**2+y**2)

        if acomp == 'z':
            amna = -(self.G*M*B*z)/(((R**2+B**2)**1.5)*(np.sqrt(z**2+H**2)))
        if acomp == 'y':
            amna = -(self.G*M*y)/((R**2+B**2)**1.5)
        if acomp == 'x':
            amna = -(self.G*M*x)/((R**2+B**2)**1.5)

        return amna

    def M31Accel(self,x,y,z,acomp):

        a_b = self.HernquistAccel(self.Mbulge,self.rbulge,x,y,z,acomp)
        a_h = self.HernquistAccel(self.Mhalo,self.rhalo,x,y,z,acomp)
        a_d = self.MiyamotoNagaiAccel(self.Mdisk,self.rd,x,y,z,acomp)

        total = a_b + a_h + a_d

        return total


    def LeapFrog(self,dt,x,y,z,vx,vy,vz):

        xhalf = x+vx*dt/2
        yhalf = y+vy*dt/2
        zhalf = z+vz*dt/2


        M31A_x = self.M31Accel(xhalf,yhalf,zhalf,'x')
        M31A_y = self.M31Accel(xhalf,yhalf,zhalf,'y')
        M31A_z = self.M31Accel(xhalf,yhalf,zhalf,'z')


        vxfull=vx+M31A_x*dt
        vyfull=vy+M31A_y*dt
        vzfull=vz+M31A_z*dt

        #calculate final position at a full time step
        xf=x+.5*(vx+vxfull)*dt
        yf=y+.5*(vy+vyfull)*dt
        zf=y+.5*(vz+vzfull)*dt


        return [xf,yf,zf,vxfull,vyfull,vzfull]

    def OrbitIntegrator(self,ti,tf,dt):

        x=self.x
        y=self.y
        z=self.z
        vx=self.vx
        vy=self.vy
        vz=self.vz


        stepper = int((tf-ti)/dt)+1
        Orbit = np.zeros((stepper,7))

        Orbit[0,1], Orbit[0,2], Orbit[0,3] = x, y, z
        Orbit[0,4], Orbit[0,5], Orbit[0,6] = vx, vy, vz
        Orbit[0,0] = ti

        #np.savetxt('M33Analytic', Orbit, header='t   x    y    z   vx   vy   vz', comments='# ',
                #fmt=['%.2f', '%.2f','%.2f','%.2f','%.2f','%.2f','%.2f'])
        t=ti+dt
        n=1
        while t <= tf:

            output = self.LeapFrog(dt,x,y,z,vx,vy,vz)


            Orbit[n,1], Orbit[n,2], Orbit[n,3] = output[0], output[1], output[2]
            Orbit[n,4], Orbit[n,5], Orbit[n,6] = output[3], output[4], output[5]

            Orbit[n,0]=t

            n = 1+n
            t = t+dt
            x,y,z = output[0], output[1], output[2]
            vx,vy,vz = output[3], output[4], output[5]

        #Save results into text file
        np.savetxt('M33Analytic', Orbit, header='t   x    y    z   vx   vy   vz', comments='# ',
                fmt=['%.2f', '%.2f','%.2f','%.2f','%.2f','%.2f','%.2f'])

        return n
