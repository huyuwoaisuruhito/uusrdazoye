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

    def get_bond_angle(self, a, o, b, f=1):
        r_NN = np.array([0, 0, -1])
        if len(set([a, o, b]))!=3:
            return 0
        OA = np.array(self.__atoms[a][1:]) - np.array(self.__atoms[o][1:])
        OB = np.array(self.__atoms[b][1:]) - np.array(self.__atoms[o][1:])
        ON = np.cross(OA, OB)
        ANGLE = np.arccos( OA.dot(OB)/(np.sqrt(OA.dot(OA)) * np.sqrt(OB.dot(OB))) )

        if f:
            if np.isnan(ANGLE):
                return 0
            elif r_NN.dot(ON)>=0:
                return ANGLE
            elif r_NN.dot(ON)<0:
                return 2*np.pi - ANGLE
            else:
                return
        else:
            if np.isnan(ANGLE):
                return 0
            elif r_NN.dot(ON)>=0:
                return ANGLE
            elif r_NN.dot(ON)<0:
                return -ANGLE
            else:
                return

    def get_dihedral_angle(self, a, b, c, d):
        r_NN = np.array([0, 0, -1])
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
        elif NN.dot(BC)>=0:
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
        self.__atoms[b] = self.__atoms[b][0:1] + list(map(lambda A, B:A + (B-A)*l/self.get_bond_length(a, b), self.__atoms[a][1:], self.__atoms[b][1:]))

    def modify_bond_angle(self, a, o, b, angle):
        if len(set([a, o, b]))!=3:
            return
        k = tuple([a, o, b])
        OA = np.array(self.__atoms[a][1:]) - np.array(self.__atoms[o][1:])
        OB = np.array(self.__atoms[b][1:]) - np.array(self.__atoms[o][1:])
        ON = self.__N_AOB.get(k, np.cross(OA, OB))
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
        else:
            raise RuntimeWarning

    def __get_unconnected_atom(self, a, b):
        unc_atom = set()

        def __(a, b, root, past):
            for i in self.__bonds[b].keys():
                if i in past: continue
                if i == a:
                    root = set()
                    return 0
                else:
                    if not __(a, i, root, past + [i]): return 0
                    root.add(i)
            return 1

        for i in self.__bonds[b].keys():
            if i == a: continue
            root = set([i])
            if __(a, i, root, [b, i]):
                unc_atom = set.union(unc_atom, root)
        
        return list(unc_atom) + [b]

    def modify_bond_length_G(self, a, b, **kw):
        if a == b:
            return
        l = kw['x']
        b_group = self.__get_unconnected_atom(a, b)
        for b_g in b_group:
            BL = self.get_bond_length(a, b)
            d_AB = list(map(lambda A, B: (B-A)*(l-BL)/BL, self.__atoms[a][1:], self.__atoms[b][1:]))
            self.__atoms[b_g] = self.__atoms[b_g][0:1] + list(map(sum, zip(self.__atoms[b_g][1:], d_AB)))

    def __get_ritation_matrix(self, v0, v1, theta):
        a, b, c = v0
        u, v, w = v1
        l = np.sqrt(sum(map(lambda x: x**2, (u, v, w))))
        u, v, w = map(lambda x: x/l, (u, v, w))
        print(v0, (u,v,w), '\ntheta:', theta)

        uu, uv, uw = u * u, u * v, u * w
        vv, vw, ww = v * v, v * w, w * w
        au, av, aw = a * u, a * v, a * w
        bu, bv, bw = b * u, b * v, b * w
        cu, cv, cw = c * u, c * v, c * w

        costheta = np.cos(theta)
        sintheta = np.sin(theta)

        m = [[0 for i in range(4)]for i in range(4)]

        m[0][0] = uu + (vv + ww) * costheta
        m[0][1] = uv * (1 - costheta) + w * sintheta
        m[0][2] = uw * (1 - costheta) - v * sintheta
        m[0][3] = 0

        m[1][0] = uv * (1 - costheta) - w * sintheta
        m[1][1] = vv + (uu + ww) * costheta
        m[1][2] = vw * (1 - costheta) + u * sintheta
        m[1][3] = 0

        m[2][0] = uw * (1 - costheta) + v * sintheta
        m[2][1] = vw * (1 - costheta) - u * sintheta
        m[2][2] = ww + (uu + vv) * costheta
        m[2][3] = 0

        m[3][0] = (a * (vv + ww) - u * (bv + cw)) * (1 - costheta) + (bw - cv) * sintheta
        m[3][1] = (b * (uu + ww) - v * (au + cw)) * (1 - costheta) + (cu - aw) * sintheta
        m[3][2] = (c * (uu + vv) - w * (au + bv)) * (1 - costheta) + (av - bu) * sintheta
        m[3][3] = 1
        
        return np.matrix(m)

    def modify_bond_angle_G(self, a, o, b, **kw):
        if len(set([a, o, b]))!=3:
            return
        delta_angle = kw['delta'] * kw['mul']
        OA = np.array(self.__atoms[a][1:]) - np.array(self.__atoms[o][1:])
        OB = np.array(self.__atoms[b][1:]) - np.array(self.__atoms[o][1:])
        ON = np.cross(OA, OB)
        M = self.__get_ritation_matrix(self.__atoms[o][1:], ON, delta_angle)

        b_group = self.__get_unconnected_atom(a, b)
        for b_g in b_group:
            pos = np.matrix(self.__atoms[b_g][1:]+[1]) * M
            self.__atoms[b_g] = self.__atoms[b_g][0:1] + list(pos.A[0][:-1])
    
    def modify_dihedral_angle_G(self, a, b, c, d, **kw):
        if len(set([a, b, c, d]))!=4:
            return
        delta_angle = kw['delta']
        while delta_angle>np.pi or delta_angle<-np.pi:
            if delta_angle>0:
                delta_angle -= 2*np.pi
            elif delta_angle<0:
                delta_angle += 2*np.pi
        delta_angle *= kw['dmul']
        BC = np.array(self.__atoms[c][1:]) - np.array(self.__atoms[b][1:])
        M = self.__get_ritation_matrix(self.__atoms[b][1:], BC, delta_angle)

        if kw['dmul']>0:
            b_group = self.__get_unconnected_atom(c, b)
        else:
            b_group = self.__get_unconnected_atom(b, c)
        for b_g in b_group:
            pos = np.matrix(self.__atoms[b_g][1:]+[1]) * M
            self.__atoms[b_g] = self.__atoms[b_g][0:1] + list(pos.A[0][:-1])

    def auto_set_O(self):
        delta = [0, 0, 0]
        sums = 0
        for s, x, y, z in self.__atoms:
            sums += PT[s]
            delta[0] += x * PT[s]
            delta[1] += y * PT[s]
            delta[2] += z * PT[s]
        if sums != 0:
            self.set_O(*[-d/sums for d in delta])
        else:
            return

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
