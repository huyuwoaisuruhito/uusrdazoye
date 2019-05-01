import copy
import numpy as np

class Molecule():

    '''分子'''

    def __init__(self):
        self.__atoms = []
        self.__bonds = []
    
    def get_atoms(self):
        return copy.deepcopy(self.__atoms)
    
    def get_bonds(self):
        return copy.deepcopy(self.__bonds)
    
    def modify_atoms(self):
        return self.__atoms
    
    def modify_bonds(self):
        return self.__bonds
    
    def get_bond_length(self, a, b):
        return np.sqrt(sum(map(lambda a, b:(a-b)**2, zip(self.__atoms[a][1:], self.__atoms[b][1:]))))
    
    def get_bond_angle(self, a, o, b):
        OA = np.array(self.__atoms[a][1:]) - np.array(self.__atoms[o][1:])
        OB = np.array(self.__atoms[b][1:]) - np.array(self.__atoms[o][1:])
        return np.arccos( OA.dot(OB)/np.sqrt(OA.dot(OA)) * np.sqrt(OB.dot(OB)) )

    def get_dihedral_angle(self, a, b, c, d):
        AB = np.array(self.__atoms[b][1:]) - np.array(self.__atoms[a][1:])
        BC = np.array(self.__atoms[c][1:]) - np.array(self.__atoms[b][1:])
        CD = np.array(self.__atoms[d][1:]) - np.array(self.__atoms[c][1:])
        N_ABC = np.cross(AB, BC)
        N_BCD = np.cross(BC, CD)
        return np.arccos( N_ABC.dot(N_BCD)/np.sqrt(N_ABC.dot(N_ABC)) * np.sqrt(N_BCD.dot(N_BCD)) )

    def modify_bond_length(self, a, b, l):
        pass
    
    def modify_bond_angle(self, a, o, b, angle):
        pass
    
    def modify_dihedral_angle(self, a, b, c, d, angle):
        pass