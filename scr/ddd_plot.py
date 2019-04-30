import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

F_RADII = {
    'H': 0.1, 'B': 0.15, 'C': 0.2, 'N': 0.2, 'O': 0.2, 'F': 0.2, 'Si': 0.3, 'P': 0.3, 'S': 0.3, 'Cl': 0.3,
}
F_COLOR = {
    'H': 'whitesmoke', 'B': 0.15, 'C': 'dimgray', 'N': 'blue', 'O': 'red', 'F': 'greenyellow', 'Si': 0.3, 'P': 'magenta', 'S': 'gold', 'Cl': 'limegreen',
}


class DDD_plot():

    '''matplotlib实现3d显示'''

    def __init__(self):
        self.fig = plt.figure(facecolor = 'mediumpurple')
    
    def init(self, molecule):
        atoms = molecule.get_atoms()
        bonds = molecule.get_bonds()
        self.ax = Axes3D(self.fig, facecolor=(0.5, 0.5, 0.797))
        self.ax.set_axis_off()

        self.plot_atoms(atoms)

        _max =max([max([abs(ll) for ll in l]) for l in zip(self.ax.get_xlim3d(), self.ax.get_ylim3d(), self.ax.get_zlim3d())])
        self.ax.set_xlim3d(-_max, _max)
        self.ax.set_ylim3d(-_max, _max)
        self.ax.set_zlim3d(-_max, _max)

        self.plot_bonds(atoms, bonds, 5/_max)
    
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
    
    def plot_bonds(self, atoms, bonding, l):
        for i, bonds in enumerate(bonding):
            A = atoms[i]
            if bonds != []:
                delta = [1, 0, 0]                
                for j, bl in bonds:
                    B = atoms[j]

                    if bl == 1:
                        theta = np.linspace(0, 1, 8)
                        x = (A[1]-B[1])*theta + B[1]
                        y = (A[2]-B[2])*theta + B[2]
                        z = (A[3]-B[3])*theta + B[3]

                        self.ax.plot(x[:5], y[:5], z[:5], color = F_COLOR[B[0]], linewidth = l)
                        self.ax.plot(x[4:], y[4:], z[4:], color = F_COLOR[A[0]], linewidth = l)
                    elif bl == 2 or bl == 1.5:    
                        if len(bonds)<2 and len(bonding[j])<2:
                            delta = [1/10/np.sqrt(l), 0, 0]
                        else:
                            if len(bonds)<2:
                                A, B = B, A
                                bonds = bonding[j]
                            n1 = np.array([A[1]-B[1], A[2]-B[2], A[3]-B[3]])
                            for j, _bl in bonds:
                                if atoms[j] != B:
                                    C = atoms[j]
                                    n2 = np.array([A[1]-C[1], A[2]-C[2], A[3]-C[3]])
                            fxl = np.cross(n2, n1)
                            if sum(fxl) == 0:
                                fxl = [0.5, 0.5, 0.5]
                            delta = np.cross(n1, fxl)
                            su = np.sqrt(delta[0]**2 + delta[1]**2 + delta[2]**2)
                            delta = [delta[i]/su/10/np.sqrt(l) for i in range(3)]

                        theta = np.linspace(0, 1, 8)
                        x = (A[1]-B[1])*theta + B[1]
                        y = (A[2]-B[2])*theta + B[2]
                        z = (A[3]-B[3])*theta + B[3]

                        if bl == 2:
                            self.ax.plot(x[:5] + delta[0], y[:5] + delta[1], z[:5] + delta[2], 
                            color = F_COLOR[B[0]], linewidth = l)
                            self.ax.plot(x[4:] + delta[0], y[4:] + delta[1], z[4:] + delta[2], 
                            color = F_COLOR[A[0]], linewidth = l)

                            self.ax.plot(x[:5] - delta[0], y[:5] - delta[1], z[:5] - delta[2], 
                            color = F_COLOR[B[0]], linewidth = l)
                            self.ax.plot(x[4:] - delta[0], y[4:] - delta[1], z[4:] - delta[2], 
                            color = F_COLOR[A[0]], linewidth = l)

                        else:
                            self.ax.plot(x[:5], y[:5], z[:5], color = F_COLOR[B[0]], linewidth = l)
                            self.ax.plot(x[4:], y[4:], z[4:], color = F_COLOR[A[0]], linewidth = l)

                            for i in range(4):
                                if i<=2:
                                    c = F_COLOR[B[0]]
                                else:
                                    c = F_COLOR[A[0]]
                                self.ax.plot(x[2*i:2*i+2] - delta[0], y[2*i:2*i+2] - delta[1], z[2*i:2*i+2] - delta[2], color = c, linewidth = l)
                    
                    elif bl == 3:
                        theta = np.linspace(0, 1, 8)
                        x = (A[1]-B[1])*theta + B[1]
                        y = (A[2]-B[2])*theta + B[2]
                        z = (A[3]-B[3])*theta + B[3]

                        n0 = np.array([0.5, 0.5, 0.5])
                        n1 = np.array([A[1]-B[1], A[2]-B[2], A[3]-B[3]])
                        delta = np.cross(n1, n0)
                        su = np.sqrt(delta[0]**2 + delta[1]**2 + delta[2]**2)
                        delta = [delta[i]/su/7.5/np.sqrt(l) for i in range(3)]

                        self.ax.plot(x[:5], y[:5], z[:5], color = F_COLOR[B[0]], linewidth = l)
                        self.ax.plot(x[4:], y[4:], z[4:], color = F_COLOR[A[0]], linewidth = l)

                        self.ax.plot(x[:5] + delta[0], y[:5] + delta[1], z[:5] + delta[2], 
                        color = F_COLOR[B[0]], linewidth = l)
                        self.ax.plot(x[4:] + delta[0], y[4:] + delta[1], z[4:] + delta[2], 
                        color = F_COLOR[A[0]], linewidth = l)

                        self.ax.plot(x[:5] - delta[0], y[:5] - delta[1], z[:5] - delta[2], 
                        color = F_COLOR[B[0]], linewidth = l)
                        self.ax.plot(x[4:] - delta[0], y[4:] - delta[1], z[4:] - delta[2], 
                        color = F_COLOR[A[0]], linewidth = l)