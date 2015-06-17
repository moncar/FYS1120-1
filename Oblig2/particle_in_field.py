from pylab import *
from mpl_toolkits.mplot3d import Axes3D

class Particle_in_field:
    def __init__(self, r0, v0, mass, charge, E_field, B_field):
        """Object particle."""
        self.r0, self.v0, self.mass, self.charge = r0, v0, mass, charge
        self.E_field, self.B_field = E_field, B_field

class Solver:
    """Class for solving with an object particle"""
    def __init__(self, particle, t0,T,dt):
        self.particle = particle;
        self.t0,self.T,self.dt = t0,T,dt

    def Cromer(self, mode=0):
        """Function for Euler-Cromer method"""
        n = int(float((self.T-self.t0)/self.dt) + 1) #number of iterations
        t = linspace(self.t0,self.T,n) #time interval grid
        
        #empty 3D array(with initial values) for position
        r = zeros((n,3)); r[0] = self.particle.r0
        v = zeros((n,3)); v[0] = self.particle.v0

        #test for which field to use
        if mode == 'E_Field':
            F = lambda dummy1,dummy2, dummy3: self.particle.charge*self.particle.E_field
        elif mode == 'B_Field':
            F = lambda vel, dummy2, dummy3: self.particle.charge*cross(vel,self.particle.B_field)
        elif mode == 'Both':
            #note: omega = qB/m
            F = lambda vel,time,test: (cos(((self.particle.charge*(norm(self.particle.B_field)))
                    /self.particle.mass)*time)
                    *self.particle.E_field*test 
                    + self.particle.charge*(cross(vel,self.particle.B_field))) 
        else:
            print "USAGE: python program.py 'options' | options: E_field, B_field, both"

        #loop through all values of t 
        for i in range(len(t)-1):
            if r[i,0] <= -0.1 or r[i,0] >= 0.1:
                """E-field is only effective for -0.1 < x > 0.1"""
                tester = 0
            else:
                tester = 1

            if mode == 'Both' and norm(r[i]-0) >= 2.6:
                """test if particle is out of cyclotrone"""
                print "Particle is out.", "t =", t[i], "v =", norm(v[i])
                break
           
            #Cromer algorithm
            a = F(v[i], t[i], tester)/self.particle.mass
            v[i+1] = v[i] + a*self.dt
            r[i+1] = r[i] + v[i+1]*self.dt

            orbit_test = False #test for finding only first orbit
            if mode == 'B_Field' and i != 0:
                """test for orbital time in B-field"""
                if (((norm(r[0] - r[i-1])) >= (norm(r[0] - r[i]))) and 
                    ((norm(r[0] - r[i])) <= (norm(r[0] - r[i+1])))):
                    
                    orbit_test == True #set test True when first orbit is found
                    print "Orbit time:", t[i]
                elif orbit_test == True:
                    """Keep loop running without printing"""
                    continue

        return t[:i], r[:i], v[:i]

    def plotter(self,t,r,v):
        """Function for plotting"""
        xkcd() #for awesomeness!

        #setup
        self.t, self.v, self.r = t,v,r
        
        subplots_adjust(right=0.77) #right-adjust subplots

        #plot position per time
        subplot(211)
        plot(t,r[:,0], t,r[:,1], t,r[:,2])
        hold('on')
        title("Position")
        ylabel("r(t)");
        legend(["x","y","z"],
                labelspacing=0.3, bbox_to_anchor=[1,1.05],
                loc=2, numpoints=1)

        #plot velocity per time
        subplot(212)
        plot(t,v[:,0], t,v[:,1], t,v[:,2])
        title("Velocity")
        ylabel("v(t)"); xlabel("t")
        legend(["vx","vy","vz"],
                labelspacing=0.3, bbox_to_anchor=[1,1.05],
                loc=2, numpoints=1)
        
        #plot position in 3d
        fig = figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(r[:,0], r[:,1], r[:,2])
        ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
        title("Position in 3D")
        show()

def _main():
    """main function"""
    #parametric values given in assigment
    r0 = (0.,0.,0.); v0_E = (0.,0.,0.); v0_B = (5.,0.,2.)     #initial values
    mass = 2.; charge = 3.;                                   #constant values for particle
    t0 = 0.; T_E = 1.; T_B = 5.; T_S = 50.; dt = 1E-04        #time values
   
    #create particle in E-field and solve using above classes
    particle_in_E_field = Particle_in_field(r0,v0_E, mass,charge, array((1.,2.,-5.)), array((0.,0.,3.)))
    solver_E = Solver(particle_in_E_field, t0,T_E,dt)
    t_E, r_E, v_E = solver_E.Cromer(mode='E_Field')
    solver_E.plotter(t_E,r_E,v_E)
     
    #plot x-position with analytical in E-field
    plot(t_E,r_E[:,0])   
    r_analytical = lambda time: r0[0] + v0_E[0]*time + (charge*5/(2*mass))*time**2
    plot(t_E,r_analytical(t_E))
    xlabel("x"); ylabel("x(t)")
    legend(["approximated","analytical"])
    show()
   
    #create particle in B-field and solve using above classes
    particle_in_B_field = Particle_in_field(r0,v0_B, mass,charge, array((1.,2.,-5.)), array((0.,0.,3.)))
    solver_B = Solver(particle_in_B_field, t0,T_B,dt)
    t_B, r_B, v_B = solver_B.Cromer(mode='B_Field')
    solver_B.plotter(t_B,r_B,v_B)
    
    #Particle in cyclotrone
    particle_in_cyclotrone = Particle_in_field(r0,(0.,0.,0.), 1.,1., array((1.,0.,0.)), array((0.,0.,1.)))
    solver_S = Solver(particle_in_cyclotrone, t0,T_S,dt)
    t_S, r_S, v_S = solver_S.Cromer(mode='Both')
    solver_S.plotter(t_S,r_S,v_S)
    
if __name__=='__main__':
    _main()
