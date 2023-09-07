from matplotlib import pyplot as plt
from matplotlib import patches as mpatches
import numpy as np 
import VTOLparam as P
# if you are having difficulty with the graphics, 
# try using one of the following backends  
# See https://matplotlib.org/stable/users/explain/backends.html
# import matplotlib
# matplotlib.use('qtagg')  # requires pyqt or pyside
# matplotlib.use('ipympl')  # requires ipympl
# matplotlib.use('gtk3agg')  # requires pyGObject and pycairo
# matplotlib.use('gtk4agg')  # requires pyGObject and pycairo
# matplotlib.use('gtk3cairo')  # requires pyGObject and pycairo
# matplotlib.use('gtk4cairo')  # requires pyGObject and pycairo
# matplotlib.use('tkagg')  # requires TkInter
# matplotlib.use('wxagg')  # requires wxPython


class VTOLanimation:
    def __init__(self):
        self.flag_init = True  # Used to indicate initialization
        # Initialize a figure and axes object
        self.fig, self.ax = plt.subplots()
        # Initializes a list of objects (patches and lines)
        self.handle = []
        # Specify the x,y axis limits
        plt.axis([-1, 20, -1, 20])
        # Draw line for the ground
        plt.plot([-10, 10], [0, 0], 'b--')
        # label axes
        plt.xlabel('z')

    def update(self, state):
        z = state[0][0]  # Horizontal position of takeoff, m
        h = state[1][0]   # Height of drone, m
        theta = state[2][0] # Angle between drone and ground, rads
        theta = (np.pi/180)*theta
        # draw plot elements: block
        plt.axis([-1, 20, -1, 20])
        self.draw_center(z, h, theta)
        self.draw_lprop(z, h, theta)
        self.draw_rprop(z, h,theta)
        self.ax.axis('equal')
        # Set initialization flag to False after first call
        if self.flag_init == True:
            self.flag_init = False

    def draw_center(self, z, h, theta):
        # specify bottom left corner of rectangle
        x = z-P.wc/2.0
        y = h-P.hc/2.0
        corner = (x, y)
        # create rectangle on first call, update on subsequent calls
        if self.flag_init is True:
            # Create the Rectangle patch and append its handle
            # to the handle list
            self.handle.append(
                mpatches.Rectangle(corner, P.wc, P.hc, rotation_point=(z,h), fc='blue', ec='black', angle=theta))
            # Add the patch to the axes
            self.ax.add_patch(self.handle[0])
        else:

            self.handle[0].set_xy(corner)  # Update patch
        
    def draw_lprop(self, z, h, theta):
        # specify centers of props
        x = [z - ((P.wc/2) + (P.d*np.cos(-theta)))]
        y = [h + (P.d*np.sin(-theta))]
        center = (x, y)
        # create circles on first call, update on previous call
        if self.flag_init is True:
            self.handle.append(
                mpatches.CirclePolygon(center, radius=0.1, resolution=15, fc='blue', ec='black')
            )
            self.ax.add_patch(self.handle[1])
        else:
            self.handle[1].xy = center

    def draw_rprop(self, z, h, theta):
        # specify centers of props
        x = [z + ((P.wc/2) + (P.d*np.cos(-theta)))]
        y = [h - (P.d*np.sin(-theta))]
        center = (x, y)
        # create circles on first call, update on previous call
        if self.flag_init is True:
            self.handle.append(
                mpatches.CirclePolygon(center, radius=0.1, resolution=15, fc='blue', ec='black')
            )
            self.ax.add_patch(self.handle[2])
        else:
            self.handle[2].xy = center
    