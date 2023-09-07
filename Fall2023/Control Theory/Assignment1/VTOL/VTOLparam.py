# Mass Spring Damper Parameter File
import numpy as np

# Physical parameters of the MSD known to the controller
mc = 1      # Mass of the drone center, kg
jc = 0.0042 # kgm^2
mr = 0.25   # Mass of right motor, kg
ml = 0.25   # Mass of left motor, kg
d = 0.3     # Distance from center to motors, m
mu = 0.1    # kg/s
g = 9.81    # gravitational accel, m/s^2

# parameters for animation
wc = 0.3     # Width of the center, m
hc = 0.3     # Height of center, m
wm = 0.3     # Width of motors, m
hm = 0.15    # Height of motors, m


# Initial Conditions
fr0 = 0.0               # ,N
fl0 = 0.0               # ,N
zv0 = 0.0               # ,m
zt0 = 0.0               # ,m
h0 = 0.0                # ,m
theta0 = 0.0            # ,rads
zvdot0 = 0.0            # ,m/s
ztdot0 = 0.0            # ,m/s
hdot0 = 0.0             # ,m/s
thetadot0 = 0.0         # ,rads

# Simulation Parameters
t_start = 0.0  # Start time of simulation
t_end = 50.0  # End time of simulation
Ts = 0.01  # sample time for simulation
t_plot = 0.1  # the plotting and animation is updated at this rate

# saturation limits
F_max = 25.0                # Max Force, N