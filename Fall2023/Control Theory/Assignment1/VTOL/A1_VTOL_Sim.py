import sys
import numpy as np
import matplotlib.pyplot as plt
import VTOLparam as P
from signalGenerator import signalGenerator
from VTOLanimation import VTOLanimation
from dataPlotter import dataPlotter
from VTOLdynamics import VTOLdynamics

# instantiate reference classes
VTOL = VTOLdynamics(alpha=0.0)
z_reference = signalGenerator(amplitude=0.5, frequency=0.1)
h_reference = signalGenerator(amplitude=0.5, frequency=0.1)
zRef = signalGenerator(amplitude=4.0, frequency=0.1, y_offset=5.0)
hRef = signalGenerator(amplitude=2.0, frequency=0.1, y_offset=2.0)
thetaRef = signalGenerator(amplitude=np.pi/8.0, frequency=0.5, y_offset=0.0)
fRef = signalGenerator(amplitude=5, frequency=0.5)
tauRef = signalGenerator(amplitude=5, frequency=0.5)


# instantiate the simulation plots and animation
#dataPlot = dataPlotter()
animation = VTOLanimation()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    while t < 7: 
        # set variables
        u = [[20], [5]]
        y = VTOL.update(u)
        # update animation
        animation.update(VTOL.state)
        #dataPlot.update(t, state, z_r, h_r, f, tau)
        print(t)
        t = t + P.t_plot
        plt.pause(0.1)

    while t >= 7 and t < 12:
        # set variables
        u = [[-30], [30]]
        y = VTOL.update(u)
        # update animation
        animation.update(VTOL.state)
        #dataPlot.update(t, state, z_r, h_r, f, tau)
        print(t)
        t = t + P.t_plot
        plt.pause(0.1)

    # set variables
    u = [[30], [20]]
    y = VTOL.update(u)
    # update animation
    animation.update(VTOL.state)
    #dataPlot.update(t, state, z_r, h_r, f, tau)
    print(t)
    t = t + P.t_plot
    plt.pause(0.1)
    
# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()