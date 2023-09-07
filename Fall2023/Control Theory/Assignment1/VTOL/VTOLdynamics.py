import sys
sys.path.append('.')
sys.path.append('/Users/benru/OneDrive/Documents/Fall 2023/Control Systems/A1 Sim/parameters')
import numpy as np 
import VTOLparam as P

class VTOLdynamics:
    def __init__(self, alpha=0.0):
        # Initial state conditions
        self.state = np.array([
            [P.zt0],  # zt initial position
            [P.h0],  # h initial value
            [P.theta0],  # theta initial value
            [P.ztdot0],  # ztdot initial velocity
            [P.hdot0],  # hdot initial
            [P.thetadot0] # theta initial radial velocity
        ])
        # simulation time step
        self.Ts = P.Ts
        # Mass of the center, kg
        self.mc = P.mc * (1.+alpha*(2.*np.random.rand()-1.))
        # , kgs^2
        self.jc = P.jc * (1.+alpha*(2.*np.random.rand()-1.))
        # Mass of right motor, kg
        self.mr = P.mr * (1.+alpha*(2.*np.random.rand()-1.))
        # Mass of left motor, kg
        self.ml = P.ml * (1.+alpha*(2.*np.random.rand()-1.))
        # Distance from center to motors, m
        self.d = P.d * (1.+alpha*(2.*np.random.rand()-1.))
        # , kg/s
        self.mu = P.mu * (1.+alpha*(2.*np.random.rand()-1.))
        # Gravitational accel, m/s^2
        self.g = P.g * (1.+alpha*(2.*np.random.rand()-1.))

        self.force_limit = P.F_max

    def update(self, u):
        # This is the external method that takes the input u at time
        # t and returns the output y at time t.
        # saturate the input force
        #u = saturate(u, self.force_limit)
        self.rk4_step(u)  # propagate the state by one time sample
        y = self.h()  # return the corresponding output
        return y
    
    def f(self, state, u):
        # Return xdot = f(x,u)
        z = state[0][0]
        h = state[1][0]
        theta = state[2][0]
        theta = (np.pi/180)*theta
        zdot = state[3][0]
        hdot = state[4][0]
        thetadot = state[5][0]
        fr = u[0][0]
        fl = u[1][0]
        # The equations of motion.
        #M = np.array([[self.mc + 2*self.mr, 0, 0],
                      #0, self.mc + 2*self.mr, 0],
                      #[0, 0, self.jc + 2*self.mr*(self.d**2)]])
        #C = np.array([[-(fr + fl)*np.sin(theta) - self.mu*zdot],
                      #[(fr + fl)*np.cos(theta) - (self.mc + 2*self.mr)*self.g],
                      #[self.d*(fr - fl)]])
        #tmp = np.linalg.inv(M) @ C
        #zddot = tmp[0][0]
        zddot = -(1/(self.mc+2*self.mr))*((fr+fl)*np.sin(theta) - self.mu*zdot)
       # hddot = tmp[1][0]
        hddot = (1/(self.mc+2*self.mr))*((fr+fl)*np.cos(theta)) - self.g
        #thetaddot = tmp[2][0]
        thetaddot = (self.d*(fr-fl))/(self.jc + 2*self.mr*(self.d**2))
        # build xdot and return
        xdot = np.array([[zdot], [hdot], [thetadot], [zddot], [hddot], [thetaddot]])
        return xdot
        

    def h(self):
        # return y = h(x)
        z = self.state[0][0]
        h = self.state[1][0]
        theta = self.state[2][0]
        y = np.array([[z], [h], [theta]])
        return y

    def rk4_step(self, u):
        # Integrate ODE using Runge-Kutta RK4 algorithm
        F1 = self.f(self.state, u)
        F2 = self.f(self.state + self.Ts / 2 * F1, u)
        F3 = self.f(self.state + self.Ts / 2 * F2, u)
        F4 = self.f(self.state + self.Ts * F3, u)
        self.state += self.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)

        
#def saturate(u, limit):
    #if abs(u) > limit:
        #u = limit*np.sign(u)
    #return u