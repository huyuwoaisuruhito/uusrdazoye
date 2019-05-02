import copy
import numpy as np

PT = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 
    'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al':13, 'Si': 14, 'P':15, 'S': 16, 
    'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 
}

class Molecule():

    '''分子'''

    def __init__(self):
        self.__charge = 0
        self.__multiplicity = 1
        self.__atoms = []
        self.__bonds = []

        self.__N_AOB = {}

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
        if a == b:
            return 0
        return np.sqrt(sum(map(lambda A, B:(A-B)**2, self.__atoms[a][1:], self.__atoms[b][1:])))

    def get_bond_angle(self, a, o, b):
        if len(set([a, o, b]))!=3:
            return 0
        OA = np.array(self.__atoms[a][1:]) - np.array(self.__atoms[o][1:])
        OB = np.array(self.__atoms[b][1:]) - np.array(self.__atoms[o][1:])
        ANGLE = np.arccos( OA.dot(OB)/(np.sqrt(OA.dot(OA)) * np.sqrt(OB.dot(OB))) )
        if np.isnan(ANGLE):
            if ANGLE<0:
                return -np.pi
            else:
                return np.pi
        else:
            return ANGLE

    def get_dihedral_angle(self, a, b, c, d):
        if len(set([a, b, c, d]))!=4:
            return 0
        AB = np.array(self.__atoms[b][1:]) - np.array(self.__atoms[a][1:])
        BC = np.array(self.__atoms[c][1:]) - np.array(self.__atoms[b][1:])
        CD = np.array(self.__atoms[d][1:]) - np.array(self.__atoms[c][1:])
        N_ABC = np.cross(AB, BC)
        N_BCD = np.cross(CD, BC)
        NN = np.cross(N_BCD, N_ABC)
        ANGLE = np.arccos( N_ABC.dot(N_BCD)/(np.sqrt(N_ABC.dot(N_ABC)) * np.sqrt(N_BCD.dot(N_BCD))) )
        if np.isnan(ANGLE):
            return 0
        elif NN.dot(BC)>0:
            return ANGLE
        elif NN.dot(BC)<0:
            return -ANGLE
        else:
            return

    def get_bond_level(self, a, b):
        return self.__bonds[a].get(b, 0)

    def modify_bond_level(self, a, b, bl):
        if a == b:
            return
        if bl != 0:
            self.__bonds[a][b] = bl
            self.__bonds[b][a] = bl
        elif bl == 0:
            try:
                del self.__bonds[a][b]
                del self.__bonds[b][a]
            except KeyError:
                pass

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
        ON = self.__N_AOB.get(k, np.cross(OA, OB))
        _ON = np.cross(OA, OB)
        if any([_ON[i]*ON[i]<0 for i in range(3)]):
            ON = -ON
        OH = np.cross(ON, OA)
        D_OA = [l/self.get_bond_length(o, a) for l in OA]
        D_OH = [l/np.sqrt(sum(map(lambda a: a**2, OH))) for l in OH]
        L_OB = self.get_bond_length(o, b)
        New_OB = list(map(sum, zip(map(lambda x: x*L_OB*np.cos(angle), D_OA), map(lambda x: x*L_OB*np.sin(angle), D_OH))))
        if not np.isnan(L_OB) and not any([np.isnan(x) for x in D_OA]) and not any([np.isnan(x) for x in D_OH]):
            self.__atoms[b] = self.__atoms[b][0:1] + list(map(sum, zip(New_OB, self.__atoms[o][1:])))
            self.__N_AOB[k] = self.__N_AOB.get(k, ON)
        else:
            raise RuntimeWarning

    def modify_dihedral_angle(self, a, b, c, d, angle):
        if len(set([a, b, c, d]))!=4:
            return
        AB = np.array(self.__atoms[b][1:]) - np.array(self.__atoms[a][1:])
        BC = np.array(self.__atoms[c][1:]) - np.array(self.__atoms[b][1:])
        CD = np.array(self.__atoms[d][1:]) - np.array(self.__atoms[c][1:])
        ANGLE = np.arccos( CD.dot(BC)/(np.sqrt(CD.dot(CD)) * np.sqrt(BC.dot(BC))) )
        N_ABC = np.cross(AB, BC)
        OH = np.cross(N_ABC, BC)
        D_N_ABC = [l/np.sqrt(sum(map(lambda a: a**2, N_ABC))) for l in N_ABC]
        D_OH = [l/np.sqrt(sum(map(lambda a: a**2, OH))) for l in OH]
        D_BC = [l/np.sqrt(sum(map(lambda a: a**2, BC))) for l in BC]
        N_P = list(map(sum, zip(map(lambda x: x*np.cos(angle), D_N_ABC), map(lambda x: x*np.sin(angle), D_OH))))
        New_CD = list(map(
            sum, zip(
                map(lambda x: x*self.get_bond_length(c, d)*np.sin(ANGLE), np.cross(D_BC, N_P)), 
                map(lambda x: x*self.get_bond_length(c, d)*np.cos(ANGLE), D_BC)
            )
        ))
        if not any([np.isnan(x) for x in New_CD]):
            self.__atoms[d] = self.__atoms[d][0:1] + list(map(sum, zip(New_CD, self.__atoms[c][1:])))
            self.__N_AOB = {}
            print(self.get_dihedral_angle(a, b, c, d))
        else:
            raise RuntimeWarning

    def modify_bond_length_G(self, a, b, l):
        pass

    def modify_bond_angle_G(self, a, o, b, angle):
        pass
    
    def modify_dihedral_angle_G(self, a, b, c, d, angle):
        pass

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
