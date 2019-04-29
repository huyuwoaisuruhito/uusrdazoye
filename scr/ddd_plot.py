import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

F_RADII = {
    'H': 0.1, 'B': 0.15, 'C': 0.2, 'N': 0.2, 'O': 0.2, 'F': 0.2, 'Si': 0.3, 'P': 0.3, 'S': 0.3, 'Cl': 0.3,
}
F_COLOR = {
    'H': 'whitesmoke', 'B': 0.15, 'C': 'dimgray', 'N': 'royablue', 'O': 'red', 'F': 'greenyellow', 'Si': 0.3, 'P': 'magenta', 'S': 'gold', 'Cl': 'limegreen',
}


class DDD_plot():

    '''matplotlib实现3d显示'''

    def __init__(self):
        self.fig = plt.figure(facecolor = 'mediumpurple')
    
    def init(self, atoms, bonding):
        self.ax = Axes3D(self.fig)
        self.ax.set_axis_off()
        self.plot_atoms(atoms)
        self.plot_bonds(atoms, bonding)
    
    def plot_atoms(self, atoms):
        for atom in atoms:
            r = F_RADII[atom[0]]
            c = F_COLOR[atom[0]]

            u = np.linspace(0, 2 * np.pi, 10)
            v = np.linspace(0, np.pi, 10)
            _x = r * np.outer(np.cos(u), np.sin(v)) + atom[1]
            _y = r * np.outer(np.sin(u), np.sin(v)) + atom[2]
            _z = r * np.outer(np.ones(np.size(u)), np.cos(v)) + atom[3]

            self.ax.plot_surface(_x, _y, _z, color=c, alpha=0.6)
    
    def plot_bonds(self, atoms, bonding):
        for i, bonds in enumerate(bonding):
            A = atoms[i]
            if bonds != []:
                for j, l in bonds:
                    B = atoms[j]
                    if l == 1:
                        theta = np.linspace(0, 1, 100)
                        x = (A[1]-B[1])*theta + B[1]
                        y = (A[2]-B[2])*theta + B[2]
                        z = (A[3]-B[3])*theta + B[3]

                        self.ax.plot(x[:51], y[:51], z[:51], color = F_COLOR[B[0]], linewidth = 5.0)
                        self.ax.plot(x[50:], y[50:], z[50:], color = F_COLOR[A[0]], linewidth = 5.0)
                        '''
                        theta = np.abs(np.pi/2 - np.arctan( (A[3]-B[3]) / np.sqrt( (A[1]-B[1])**2 + (A[2]-B[2])**2 ) ))#极角
                        phi = np.arctan( (A[2]-B[2]) / (A[1]-B[1]) )#方位角

                        u = np.linspace(0, 2*np.pi, 20)
                        h = np.linspace(0, 1, 100)
                        
                        x = np.outer(np.sin(u), np.ones(len(h)))
                        y = np.outer(np.cos(u), np.ones(len(h)))
                        _z = np.outer(np.ones(len(u)), h) + x * np.sin(theta)
                        x = x * np.cos(theta)

                        _x, _y = x*np.cos(phi) - y*np.sin(phi), x*np.sin(phi) + y*np.cos(phi)
                        l = np.sqrt( (A[1]-B[1])**2 + (A[2]-B[2])**2 + (A[3]-B[3])**2)
                        _x, _y, _z = 0.05*_x, 0.05*_y, l*_z

                        self.ax.plot_surface(_x, _y, _z, color='r')
                        '''