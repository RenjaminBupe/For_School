
import sys
sys.path.append('.')# one directory up
import numpy as np
import time 
from math import cos, sin, tan
import scipy.linalg as linalg
import parameters.simulation_parameters as SIM
import matplotlib.pyplot as plt
from viewers.UAV_animation import UAV_animation
from tools.signalGenerator import signalGenerator
from tools.sliders import sliders
from dynamics.UAVDynamics_Assignment_2 import UAVdynamics

UAV = UAVdynamics(alpha=0.0)
state=np.array([[0.], [0.], [-1.], [0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.]])
UAV_anim=UAV_animation(state, scale=5)
#my_slider=sliders()
   
fref = signalGenerator(1000., 2.)
mref = signalGenerator(30., 2)

# initialize the simulation time
sim_time = SIM.start_time

# initialize plots


# initialize plot lists
pn_p = [0]
pe_p = [0]
pd_p = [0]


# main simulation loop
print("Press Command-Q to exit...")
pn=state[0][0]
pe=state[1][0]
pd=state[2][0]
u = state[3][0]
v = state[4][0]
w = state[5][0]
phi=state[6][0]
theta=state[7][0]
psi=state[8][0]
p = state[9][0]
q = state[10][0]
r = state[11][0]

while sim_time < SIM.end_time:
    if sim_time <= 0.5:
        fx = 0
        fy = fref.sin(sim_time)
        fz = 0
        l = 0
        m = 0
        n = 0
    elif sim_time <= 1.0:
        fy = 0
        fx = fref.sin(sim_time)
    elif sim_time <= 1.5:
        fx = 0
        fz = -fref.sin(sim_time)
    elif sim_time <= 2.0:
        fz = 0
        n = mref.sin(sim_time)
    elif sim_time <= 2.5:
        n = 0
        m = mref.sin(sim_time)
    elif sim_time <= 3.0:
        m = 0
        l = mref.sin(sim_time)
    else:
        l = 0

    
    y = UAV.update(fx, fy, fz, l, m, n)
    pn, pe, pd, u, v, w, phi, theta, psi, p, q, r = y.flatten()
    UAV_anim.update(pn, pe, pd, phi, theta, psi) # -pd for height
    

    # -------increment time-------------d
    plt.pause(0.1)
    sim_time += SIM.ts_simulation

    



