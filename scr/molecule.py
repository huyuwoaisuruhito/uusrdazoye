import copy
import numpy as np

class Molecule():

    '''分子'''

    def __init__(self):
        self.__charge = 0
        self.__multiplicity = 1
        self.__atoms = []
        self.__bonds = []

        self.__N_A = {}

    def copy(self):
        return copy.deepcopy(self)

    def set_charge(self, c):
        self.__charge = c
    
    def set_multiplicity(self, m):
        self.__multiplicity = m
    
    def get_charge(self, c):
        return self.__charge
    
    def get_multiplicity(self, m):
        return self.__multiplicity

    def get_atoms(self):
        return copy.deepcopy(self.__atoms)

    def get_bonds(self):
        return copy.deepcopy(self.__bonds)

    def modify_atoms(self):
        return self.__atoms

    def modify_bonds(self):
        return self.__bonds

    def get_bond_length(self, a, b):
        return np.sqrt(sum(map(lambda A, B:(A-B)**2, self.__atoms[a][1:], self.__atoms[b][1:])))

    def get_bond_angle(self, a, o, b):
        OA = np.array(self.__atoms[a][1:]) - np.array(self.__atoms[o][1:])
        OB = np.array(self.__atoms[b][1:]) - np.array(self.__atoms[o][1:])
        return np.arccos( OA.dot(OB)/(np.sqrt(OA.dot(OA)) * np.sqrt(OB.dot(OB))) )

    def get_dihedral_angle(self, a, b, c, d):
        AB = np.array(self.__atoms[b][1:]) - np.array(self.__atoms[a][1:])
        BC = np.array(self.__atoms[c][1:]) - np.array(self.__atoms[b][1:])
        CD = np.array(self.__atoms[d][1:]) - np.array(self.__atoms[c][1:])
        N_ABC = np.cross(AB, BC)
        N_BCD = np.cross(BC, CD)
        return np.arccos( N_ABC.dot(N_BCD)/(np.sqrt(N_ABC.dot(N_ABC)) * np.sqrt(N_BCD.dot(N_BCD))) )

    def modify_bond_length(self, a, b, l):
        if a == b:
            return
        BL = self.get_bond_length(a, b)
        self.__atoms[b] = self.__atoms[b][0:1] + list(map(lambda A, B:A + (B-A)*l/self.get_bond_length(a, b), self.__atoms[a][1:], self.__atoms[b][1:]))

    def modify_bond_angle(self, a, o, b, angle):
        if len(set([a, o, b]))!=3:
            return
        OA = np.array(self.__atoms[a][1:]) - np.array(self.__atoms[o][1:])
        OB = np.array(self.__atoms[b][1:]) - np.array(self.__atoms[o][1:])
        k = tuple(sorted([a, b, o]))
        OH = self.__N_A.get(k, np.cross(np.cross(OA, OB), OA))
        if OH[0]<0:
            OH = -OH
        D_OA = [l/self.get_bond_length(o, a) for l in OA]
        D_OH = [l/np.sqrt(sum(map(lambda a: a**2, OH))) for l in OH]
        L_OB = self.get_bond_length(o, b)
        New_OB = list(map(sum, zip(map(lambda x: x*L_OB*np.cos(angle), D_OA), map(lambda x: x*L_OB*np.sin(angle), D_OH))))
        if not np.isnan(L_OB) and not any([np.isnan(x) for x in D_OA]) and not any([np.isnan(x) for x in D_OH]):
            self.__atoms[b] = self.__atoms[b][0:1] + list(map(sum, zip(New_OB, self.__atoms[o][1:])))
            self.__N_A[k] = OH
        else:
            raise RuntimeWarning

    def modify_dihedral_angle(self, a, b, c, d, angle):
        if len(set([a, b, c, d]))!=4:
            return
        AB = np.array(self.__atoms[b][1:]) - np.array(self.__atoms[a][1:])
        BC = np.array(self.__atoms[c][1:]) - np.array(self.__atoms[b][1:])
        CD = np.array(self.__atoms[d][1:]) - np.array(self.__atoms[c][1:])
        N_ABC = np.cross(AB, BC)
        OH = np.cross(BC, N_ABC)
        D_N_ABC = [l/np.sqrt(sum(map(lambda a: a**2, N_ABC))) for l in N_ABC]
        D_OH = [l/np.sqrt(sum(map(lambda a: a**2, OH))) for l in OH]
        D_BC = [-l/np.sqrt(sum(map(lambda a: a**2, BC))) for l in BC]
        N_P = list(map(sum, zip(map(lambda x: x*np.cos(angle), D_N_ABC), map(lambda x: x*np.sin(angle), D_OH))))
        New_CD = list(map(
            sum, zip(
                map(lambda x: x*self.get_bond_length(c, d)*np.sin(self.get_bond_angle(b, c, d)), np.cross(D_BC, N_P)), 
                map(lambda x: x*self.get_bond_length(c, d)*np.cos(self.get_bond_angle(b, c, d)), D_BC)
            )
        ))
        if not any([np.isnan(x) for x in New_CD]):
            self.__atoms[d] = self.__atoms[d][0:1] + list(map(sum, zip(New_CD, self.__atoms[c][1:])))
            self.__N_A = {}
        else:
            raise RuntimeWarning

    def set_O(self, x, y, z):
        self.__atoms = [ [atom[0], atom[1]+x, atom[2]+y, atom[3]+z] for atom in self.__atoms ]

    def add_atom(self, a):
        pass

    def add_bond(self, a, b):
        pass

    def test(self):
        for i in range(360):
            self.modify_bond_angle(0, 1, 2, np.pi*i/18)
            #self.modify_dihedral_angle(2, 0, 1, 3, np.pi*i/18 - np.pi)
            print(self.__atoms[3][1:])
