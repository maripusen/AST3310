import FVis3 as fvis
import numpy.ma as ma
import numpy as np


"""
THIS IS THE ON AND OFF "BUTTON" FOR THE GAUSSION PERTURBATION 1= ON, 0=OFF
"""
perturb = 1
#a secondary term for the introduction of two other gaussian perturbations, if viewing several perturbations remember to set viewing time to 120 seconds
multiperturb =0

#A few basic constants
g = -(6.67408 * 10**(-11) * 1.989*10**30) /(6.96*10**8)**2
gamma = 5/3
mu = 0.61
m_u = 1.66*10**(-27)
k = 1.38065*10**(-23)
kb  = 5.670374419 * 10**(-8)
class convection():

    def __init__(self):
        """
        define variables
        """
        self.nx = 300
        self.ny = 100
        self.X = 12*10**6
        self.Y = 4*10**6
        self.dx = self.X/self.nx
        self.dy = self.Y/self.ny
        self.dt = 0.01
        self.u = np.zeros((self.nx,self.ny))
        self.w = np.zeros((self.nx,self.ny))
        self.rho = np.zeros((self.nx,self.ny))
        self.e = np.zeros((self.nx,self.ny))
        self.P = np.zeros((self.nx,self.ny))
        self.T = np.zeros((self.nx,self.ny))

    def initialise(self):
        """
        initialise temperature, pressure, density and internal energy
        """
        self.T[:,0] = 5770
        self.P[:,0] = 1.8*10**8
        self.e[:,0] = self.P[:,0] /(gamma - 1)
        self.rho[:,0] = self.P[:,0] *mu * m_u /(k*self.T[:,0])
        self.T[:] = 5770 -(mu * m_u * g *(2/5 + 0.001)*np.linspace(0,self.Y,int(self.Y/self.dy)))/k
        #Introducing one or more gaussian perturbations
        if perturb == 1:
            A       = 10000
            x       = np.linspace(0, self.nx - 1, self.nx)
            x_      = (self.nx - 1) / 2
            sig_x   = 10
            y_      = (self.ny - 1) / 2
            sig_y   = 10
            gauss = np.zeros((self.nx,self.ny))
            gauss2 = np.zeros((self.nx,self.ny))
            gauss3 = np.zeros((self.nx,self.ny))
            for j in range(self.ny):
                gauss[:, j] = A * np.exp( - ( (x - x_) ** 2 / (2 * sig_x ** 2) + (j - y_) ** 2 / (2 * sig_y ** 2) ) )
                if multiperturb == 1:
                    gauss2[:,j] = A * np.exp( - ( (x - x_/2) ** 2 / (2 * sig_x ** 2) + (j - y_) ** 2 / (2 * sig_y ** 2) ) )
                    gauss3[:,j] = A * np.exp( - ( (x - 2*x_) ** 2 / (2 * sig_x ** 2) + (j - y_) ** 2 / (2 * sig_y ** 2) ) )
                else:
                    gauss2 = 0
                    gauss3 = 0
            self.T[:] =self.T + gauss + gauss2 + gauss3
        self.P[:] = 1.8*10**8*(self.T[:,:]/5770)**(1/0.401)
        self.rho[:] = self.P[:,:] *mu * m_u /(k*self.T[:,:])
        self.e[:] = self.P[:,:]/(gamma -1)

    def timestep(self):
        """
        calculate timestep
        """
        #defining a change parameter
        p = 0.1

        #Extracting the flow velocities
        ux, uy = np.array(np.where(self.u!=0))
        wx, wy = np.array(np.where(self.w!=0))

        #Calculating the relative change in the flow velocities
        if len(ux) == 0:
            reu = 0
        else:
            reu = np.abs((self.drhou[ux[:],uy[:]]/self.rho[ux[:],uy[:]])/self.u[ux[:],uy[:]])
        if len(wx) == 0:
            rew = 0
        else:
            rew = np.abs((self.drhow[wx[:],wy[:]]/self.rho[wx[:],wy[:]])/self.u[wx[:],wy[:]])

        #Calculating the relative change for the other variables
        rerho = np.abs(self.drho/self.rho)
        ree = np.abs(self.de/self.e)
        rex = np.abs(self.u/self.dx)
        rey = np.abs(self.w/self.dy)

        #Extracting the maximas
        maxrho = np.max(rerho)
        maxe = np.max(ree)
        maxu = np.max(reu)
        maxw = np.max(rew)
        maxx = np.max(rex)
        maxy = np.max(rey)

        Maximilliano = [maxrho,maxe,maxu,maxw,maxx,maxy]
        #Maximilliano = [maxrho,maxe,maxx,maxy]
        #Extracting THE maximum
        dt = max(Maximilliano)
        #Then setting the actual timestep with some boundaries and conditions to make sure they dont get too large or small
        if dt == 0:
            self.dt == 0.01
        else:
            self.dt = p/dt
            if self.dt < 0.01:
                self.dt = 0.01
            elif self.dt > 10:
                self.dt = 10

    def boundary_conditions(self):
        """
        boundary conditions for energy, density and velocity
        """
        self.w[:,0] = 0
        self.w[:,-1] = 0

        self.u[:,0] = (4*self.u[:,1]- self.u[:,2])/3
        self.u[:,-1] = (4*self.u[:,-2]- self.u[:,-3])/3

        self.e[:,0] = (- self.e[:,2] + 4 *self.e[:,1])/(3-2 * mu *m_u* g*self.dy /(k*self.T[:,0]))
        self.e[:,-1] = (- self.e[:,-3] + 4 * self.e [:,-2])/(3+2 * mu*m_u*g*self.dy/(k * self.T[:,-1]))

        self.rho[:,0] = self.e[:,0] * 2/3 *mu * m_u /(k*self.T[:,0])
        self.rho[:,-1] = self.e[:,-1] * 2/3 *mu *m_u/(k* self.T[:,-1])

    def central_x(self,func):
        """
        central difference scheme in x-direction
        """
        forward = np.roll(func, 1, axis=0)
        backward = np.roll(func, -1, axis=0)
        return (forward - backward)/(2*self.dx)

    def central_y(self,func):
        """
        central difference scheme in y-direction
        """
        forward = np.roll(func, 1, axis=1)
        backward = np.roll(func, -1, axis=1)

        return (forward - backward)/(2*self.dy)

    def upwind_x(self,func,u):
        """
        upwind difference scheme in x-direction
        """

        funret = np.zeros((self.nx,self.ny))
        forward = np.roll(func, 1, axis=0)
        backward = np.roll(func, -1, axis=0)

        pos = np.where(u >= 0)
        neg = np.where(u < 0)


        funret[pos[0],pos[1]] = (func[pos[0],pos[1]]-backward[pos[0],pos[1]])/self.dx
        funret[neg[0],neg[1]] = (forward[neg[0],neg[1]]-func[neg[0],neg[1]])/self.dx
        return funret

    def upwind_y(self,func,u):
        """
        upwind difference scheme in y-direction
        """
        funret = np.zeros((self.nx,self.ny))
        forward = np.roll(func, 1, axis=1)
        backward = np.roll(func, -1, axis=1)

        pos = np.where(u >= 0)
        neg = np.where(u < 0)


        funret[pos[0],pos[1]] = (func[pos[0],pos[1]]-backward[pos[0],pos[1]])/self.dy
        funret[neg[0],neg[1]] = (forward[neg[0],neg[1]]-func[neg[0],neg[1]])/self.dy
        return funret

    def hydro_solver(self):
        """
        hydrodynamic equations solver
        """
        prevrho = self.rho
        prevu = self.u
        prevw = self.w
        preve = self.e
        prevp = self.P


        #These are not entirely correct, they must be discrete, they are now very discrete
        self.drho = -prevrho * (self.central_x(prevu)+ self.central_y(prevw))- prevu* self.upwind_x(prevrho,prevu)-prevw * self.upwind_y(prevrho,prevw)

        self.drhou = prevrho*prevu*(self.upwind_x(prevu,prevu)+ self.upwind_y(prevw,prevu))- prevu*self.upwind_x(prevrho*prevu,prevu) - prevw* self.upwind_y(prevrho*prevu,prevw)  -self.central_x(prevp)
        self.drhow = prevrho*prevw*(self.upwind_y(prevw,prevw)+ self.upwind_x(prevu,prevw))- prevw*self.upwind_y(prevrho*prevw,prevw) - prevu* self.upwind_x(prevrho*prevw,prevu)  -self.central_y(prevp) + prevrho *g
        self.de = -preve * (self.central_x(prevu)+ self.central_y(prevw))- prevu * self.upwind_x(preve,prevu)- prevw * self.upwind_y(preve,prevw)- prevp*(self.central_x(prevu)+ self.central_y(prevw))

        #Update main variables and the time step
        self.timestep()
        self.rho[:] = prevrho + self.drho * self.dt
        self.u[:] = (prevrho * prevu + self.drhou*self.dt)/(self.rho)
        self.w[:] = (prevrho * prevw + self.drhow*self.dt)/(self.rho)
        self.e[:] = preve + self.de * self.dt

        #Impose boundaries
        self.boundary_conditions()
        self.P[:] = preve*(2/3)
        self.T[:] = (prevp*mu * m_u/(k*self.rho))
        return self.dt

a = convection()
init1 = a.initialise()

vis = fvis.FluidVisualiser()

#vis.save_data(200,a.hydro_solver, rho=np.rot90(a.rho,1), u=np.rot90(a.u,1),w = np.rot90(a.w,1), P=np.rot90(a.P,1),T = np.rot90(a.T,1),e = np.rot90(a.e,1), sim_fps=1.0)
#vis.animate_2D("T",folder = "FVis_output_2020-05-18_20-01")
vis.animate_2D("rho",folder = "singularperturbation")

if __name__ == "__main__":
    pi = 3
