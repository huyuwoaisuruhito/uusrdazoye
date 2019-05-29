"""
分子力学能量优化计算

算法：梯度下降法
力场：amber03
"""


import numpy as np
from mm.bonding_information import Bondingmaps
from mm.atomtype import Atoms
from scr.molecule import Molecule


class DataTypeError(TypeError):
    pass

class Computing:
    def __init__(self, mainwindow, molecule):
        self._bondingmap = Bondingmaps(molecule)
        
        self._atomtypes = Computing.__gettypes(self._bondingmap)
        self.sites = [np.array(list(map(lambda lenth: lenth, i[1:]))) for i in self._bondingmap.sites]
        self._newsites = self.sites[:]
        
        self._parent = mainwindow
        self._molecule = molecule
        
        self._bonds = []
        self._bonddatas = dict()
        self._angles = []
        self._angledatas = dict()
        self._dihedrals = []
        self._dihedraldatas = dict()
        self._nonbonds = []
        self._nonbonddatas = dict()
        self.__getdatas()
        
        self._forces = []
        self._step = 0.01
        self._potential = 0
        self._oldpotential = 0
        self.message = ""
        self._times = 0
        
    def run(self):
        self.__calc()
        self._oldpotential = self._potential
        self.__printstart()
        self._parent.event_generate("<<Changelog>>")
        self._parent.event_generate("<<Changeview>>")
        
        while max([np.sqrt(force.dot(force)) for force in self._forces]) > 0.00005:
            
            self.__calc()
            if self._times % 5 == 1:
                atomsites = self._molecule.modify_atoms()
                for i, site in enumerate(self.sites):
                    atomsites[i][1:] = list(site)
                
                self.__printstep()
                self._parent.event_generate("<<Changelog>>")
                self._parent.event_generate("<<Changeview>>")
                self._oldpotential = self._potential
        
        self.__printfinish()
        self._parent.event_generate("<<Changelog>>")
        self._parent.event_generate("<<Changeview>>")
    
    def __printstart(self):
        self.message = "开始计算\r\n\r\n---------------------------------------"
    
    def __printstep(self):
        stepdatas = (self._times//5, self._potential, self._oldpotential - self._potential)
        self.message = "\n计算轮次: {0:d}\n总能量:\t{1:.4f} kJ/mol\n能量下降:\t{2:.4f} kJ/mol".format(*stepdatas)
    
    def __printfinish(self):
        self.message = "\n总能量:\t{0:.4f} kJ/mol\n共计算 {1:d} 次 \n\r\n计算结束\r\n\r\n---------------------------------------".format(self._potential, self._times)
    
    def __calc(self):
        #计算新坐标对应的势能
        def potential_bond(array1, array2, bonddata):
            arrayr = array1 - array2
            return bonddata[2] * (np.sqrt(arrayr.dot(arrayr)) - bonddata[1])**2 / 2
        
        def potential_angle(array1, array2, array3, angledata):
            aleft = array1 - array2
            aright = array3 - array2
            theta = np.arccos(aleft.dot(aright)/np.sqrt(aleft.dot(aleft) * aright.dot(aright)))
            return angledata[1] * (theta - angledata[0])**2 / 2
        
        def potential_dihedral(array1, array2, array3, array4, dihedraldata):
            a1 = np.cross((array2 - array1), (array3 - array2))
            a2 = np.cross((array3 - array2), (array4 - array3))
            phi = np.arccos(a1.dot(a2)/np.sqrt(a1.dot(a1) * a2.dot(a2)))
            return np.sum([data[2]*np.cos(data[0]*phi - data[1]) for data in dihedraldata[1]])
        
        def potential_nonbond(array1, array2, nonbonddata):
            arrayr = array1 - array2
            rd = nonbonddata[0] / np.sqrt(arrayr.dot(arrayr))
            return nonbonddata[1]*(rd**12 - 2*rd**6)
        
        newpotential = np.sum([potential_bond(self._newsites[b[0]], self._newsites[b[1]], self._bonddatas[b]) for b in self._bonds])
        newpotential += np.sum([potential_angle(self._newsites[ag[0]], self._newsites[ag[1]], self._newsites[ag[2]], self._angledatas[ag])
                           for ag in self._angles])
        newpotential += np.sum([potential_dihedral(*[self._newsites[i] for i in dh], self._dihedraldatas[dh])
                           for dh in self._dihedrals])
        newpotential += np.sum([potential_nonbond(self._newsites[nb[0]], self._newsites[nb[1]], self._nonbonddatas[nb])
                            for nb in self._nonbonds])
        
        if self._times == 0:
            self._potential = newpotential
        else:
            if newpotential < self._potential:
                self._step *= 1.2
                self._potential = newpotential
                self.sites = self._newsites[:]
            else:
                self._step *= 0.6
        
        #计算力（势能偏导）以得到一组新坐标
        def fbond(num1, array1, array2, bonddata):
            arrayr = array2 - array1
            r = np.sqrt(arrayr.dot(arrayr))
            return (bonddata[2] * (r - bonddata[1]) / r) * arrayr
        
        def fangle_c(num2, array1, array2, array3, angledata):
            dpotential_x = (potential_angle(array1, array2+np.array((0.00001, 0, 0)), array3, angledata)
                            - potential_angle(array1, array2, array3, angledata))
            dpotential_y = (potential_angle(array1, array2+np.array((0, 0.00001, 0)), array3, angledata)
                            - potential_angle(array1, array2, array3, angledata))
            dpotential_z = (potential_angle(array1, array2+np.array((0, 0, 0.00001)), array3, angledata)
                            - potential_angle(array1, array2, array3, angledata))
            return - np.array((dpotential_x, dpotential_y, dpotential_z)) / 0.00001
        
        def fangle_s(num1, array1, array2, array3, angledata):
            dpotential_x = (potential_angle(array1+np.array((0.00001, 0, 0)), array2, array3, angledata)
                            - potential_angle(array1, array2, array3, angledata))
            dpotential_y = (potential_angle(array1+np.array((0, 0.00001, 0)), array2, array3, angledata)
                            - potential_angle(array1, array2, array3, angledata))
            dpotential_z = (potential_angle(array1+np.array((0, 0, 0.00001)), array2, array3, angledata)
                            - potential_angle(array1, array2, array3, angledata))
            return - np.array((dpotential_x, dpotential_y, dpotential_z)) / 0.00001
        
        def fdihedral_c(num2, array1, array2, array3, array4, dihedraldata):
            dpotential_x = ((potential_dihedral(array1, array2+np.array((0.00001, 0, 0)), array3, array4, dihedraldata))
                            -(potential_dihedral(array1, array2, array3, array4, dihedraldata)))
            dpotential_y = ((potential_dihedral(array1, array2+np.array((0, 0.00001, 0)), array3, array4, dihedraldata))
                            -(potential_dihedral(array1, array2, array3, array4, dihedraldata)))
            dpotential_z = ((potential_dihedral(array1, array2+np.array((0, 0, 0.00001)), array3, array4, dihedraldata))
                            -(potential_dihedral(array1, array2, array3, array4, dihedraldata)))
            return - np.array((dpotential_x, dpotential_y, dpotential_z)) / 0.00001
        
        def fdihedral_s(num1, array1, array2, array3, array4, dihedraldata):
            dpotential_x = ((potential_dihedral(array1+np.array((0.00001, 0, 0)), array2, array3, array4, dihedraldata))
                            -(potential_dihedral(array1, array2, array3, array4, dihedraldata)))
            dpotential_y = ((potential_dihedral(array1+np.array((0, 0.00001, 0)), array2, array3, array4, dihedraldata))
                            -(potential_dihedral(array1, array2, array3, array4, dihedraldata)))
            dpotential_z = ((potential_dihedral(array1+np.array((0, 0, 0.00001)), array2, array3, array4, dihedraldata))
                            -(potential_dihedral(array1, array2, array3, array4, dihedraldata)))
            return - np.array((dpotential_x, dpotential_y, dpotential_z)) / 0.00001
        
        def fnonbond(num1, array1, array2, nonbonddata):
            arrayr = array2 - array1
            r = np.sqrt(arrayr.dot(arrayr))
            rd = nonbonddata[0] / r
            return (12 * nonbonddata[1]*(rd**12 - rd**6) / r**2) * arrayr
        
        forces = np.array([(0.0, 0.0, 0.0) for i in range(len(self.sites))])
        for bond in self._bonds:
            forces[bond[0]] += fbond(bond[0], self.sites[bond[0]], self.sites[bond[1]], self._bonddatas[bond])
            forces[bond[1]] += fbond(bond[1], self.sites[bond[1]], self.sites[bond[0]], self._bonddatas[bond])
        
        for angle in self._angles:
            forces[angle[0]] += fangle_s(angle[0], *[self.sites[at] for at in angle], self._angledatas[angle])
            forces[angle[2]] += fangle_s(angle[2], *[self.sites[at] for at in reversed(angle)], self._angledatas[angle])
            forces[angle[1]] += fangle_c(angle[1], *[self.sites[at] for at in angle], self._angledatas[angle])
        
        for dh in self._dihedrals:
            forces[dh[0]] += fdihedral_s(dh[0], *[self.sites[at] for at in dh], self._dihedraldatas[dh])
            forces[dh[3]] += fdihedral_s(dh[3], *[self.sites[at] for at in reversed(dh)], self._dihedraldatas[dh])
            forces[dh[1]] += fdihedral_c(dh[1], *[self.sites[at] for at in dh], self._dihedraldatas[dh])
            forces[dh[2]] += fdihedral_c(dh[2], *[self.sites[at] for at in reversed(dh)], self._dihedraldatas[dh])
        
        for nonbond in self._nonbonds:
            forces[nonbond[0]] += fnonbond(nonbond[0], *[self.sites[at] for at in nonbond], self._nonbonddatas[nonbond])
            forces[nonbond[1]] += fnonbond(nonbond[1], *[self.sites[at] for at in reversed(nonbond)], self._nonbonddatas[nonbond])
        
        maxforce = max([np.sqrt(force.dot(force)) for force in forces])
        if maxforce >= 0.00077:
            self._forces = forces * (self._step / maxforce)
        else:
            self._forces = forces

        self._newsites += self._forces
        self._times += 1
    
    
    def __getdatas(self):
        
        ori_bond_datas = [line.split() for line in open("datas\\bond_data.txt", "r")]
        ori_angle_datas = [line.split() for line in open("datas\\angle_data.txt", "r")]
        ori_dihedral_datas = [line.split() for line in open("datas\\dihedral_data.txt", "r")]
        ori_vdw_datas = [line.split() for line in open("datas\\vdw_data.txt", "r")]
        
        #列出键列表以及键长计算参数
        for a, bonding in enumerate(self._bondingmap.bondinglist):
            for bond in bonding:
                self._bonds.append((a, bond[0]))
                key = (self._atomtypes[a], self._atomtypes[bond[0]])
                
                for ori_data in ori_bond_datas:
                    if (key[0] == ori_data[0] and key[1] == ori_data[1]) or (key[0] == ori_data[1] and key[1] == ori_data[0]):
                        self._bonddatas[(a, bond[0])] = (bond[1], float(ori_data[3])*10, float(ori_data[4])/100)
                        break
                else:
                    raise DataTypeError("can't find data for this bond: {0}".format(str((a, bond[0]))))
        
        #列出角列表以及相关参数
        for cnum, center in enumerate(self._atomtypes):
            if len(self._bondingmap.relations[cnum]) > 1:
                for lnum, left in enumerate(self._bondingmap.relations[cnum]):
                    for right in self._bondingmap.relations[cnum][lnum+1:]:
                        self._angles.append((left, cnum, right))
                        key = (self._atomtypes[left], center, self._atomtypes[right])
                        
                        for ori_data in ori_angle_datas:
                            if ((key[0] == ori_data[0] and key[1] == ori_data[1] and key[2] == ori_data[2]) or
                                (key[0] == ori_data[2] and key[1] == ori_data[1] and key[2] == ori_data[0])
                            ):
                                self._angledatas[(left, cnum, right)] = (float(ori_data[4])*(np.pi)/180, float(ori_data[5]))
                                break
                            if ((key[0] == ori_data[0] and key[1] == ori_data[1] and key[2][0] == ori_data[2][0] == "H") or
                                (key[0] == ori_data[2] and key[1] == ori_data[1] and key[2][0] == ori_data[0][0] == "H") or
                                (key[2] == ori_data[0] and key[1] == ori_data[1] and key[0][0] == ori_data[2][0] == "H")
                            ):
                                self._angledatas[(left, cnum, right)] = (float(ori_data[4])*(np.pi)/180, float(ori_data[5]))
                                break
                        else:
                            raise DataTypeError("can't find data for this angle: {0}".format(str((left, cnum, right))))
        
        #列出二面角列表以及相关参数
        def append_dihedral(dihedral, ori_data):
            nonlocal find
            find += 1
            if ori_data[4] == "4":
                self._dihedraldatas[dihedral] = (4, [(int(ori_data[7]), float(ori_data[5])*(np.pi)/180, float(ori_data[6]))])
            if ori_data[4] == "9":
                if dihedral not in self._dihedraldatas:
                    self._dihedraldatas[dihedral] = (9, [])
                self._dihedraldatas[dihedral][1].append((int(ori_data[7]), float(ori_data[5])*(np.pi)/180, float(ori_data[6])))
        
        for cnum1, bonding in enumerate(self._bondingmap.bondinglist):
            for cbond in bonding:
                cnum2 = cbond[0]
                for left in self._bondingmap.relations[cnum1]:
                    if left == cnum2:
                        continue
                    for right in self._bondingmap.relations[cnum2]:
                        if right == cnum1:
                            continue
                        
                        dihedral = (left, cnum1, cnum2, right)
                        self._dihedrals.append(dihedral)
                        key = tuple(self._atomtypes[i] for i in (dihedral))
                        
                        find = 0
                        for ori_data in ori_dihedral_datas:
                            if ((key[0] == ori_data[0] and key[1] == ori_data[1] and key[2] == ori_data[2] and key[3] == ori_data[3]) or
                                (key[3] == ori_data[0] and key[2] == ori_data[1] and key[1] == ori_data[2] and key[0] == ori_data[3])
                            ):
                                append_dihedral(dihedral, ori_data)
                                continue
                            if find:
                                break
                            
                            if ori_data[0] == "X":
                                if ((key[1] == ori_data[1] and key[2] == ori_data[2] and key[3] == ori_data[3]) or
                                    (key[2] == ori_data[1] and key[1] == ori_data[2] and key[0] == ori_data[3])
                                ):
                                    append_dihedral(dihedral, ori_data)
                                    continue
                                if find:
                                    break
                                
                                if ori_data[3] == "X":
                                    if ((key[1] == ori_data[1] and key[2] == ori_data[2]) or
                                        (key[1] == ori_data[2] and key[2] == ori_data[1])
                                    ):
                                        append_dihedral(dihedral, ori_data)
                                        continue
                                if find:
                                    break
                                        
                                if ori_data[1] == "X":
                                    if ((key[2] == ori_data[2] and key[3] == ori_data[3]) or
                                        (key[1] == ori_data[2] and key[0] == ori_data[3])
                                    ):
                                        append_dihedral(dihedral, ori_data)
                                        continue
                                if find:
                                    break
                        
                        else:
                            raise DataTypeError("can't find data for this dihedral: {0}".format(str(dihedral)))
        
        #非键作用
        for numa, typea in enumerate(self._atomtypes):
            for ori_data_a in ori_vdw_datas:
                if typea != ori_data_a[0]:
                    continue
                ra, epa = float(ori_data_a[5]), float(ori_data_a[6])
                for numc, typeb in enumerate(self._atomtypes[numa+1:]):
                    numb = numc + numa + 1
                    if ((numa, numb) in self._bonds or (numb, numa) in self._bonds):
                        continue
                    relations_13 = [(i[0], i[2]) for i in self._angles]
                    relations_14 = [(i[0], i[3]) for i in self._dihedrals]
                    if (((numa, numb) in relations_13 or (numb, numa) in relations_13) or
                        ((numa, numb) in relations_14 or (numb, numa) in relations_14)
                    ):
                        continue
                    for ori_data_b in ori_vdw_datas:
                        if typeb != ori_data_b[0]:
                            continue
                        rb, epb = float(ori_data_b[5]), float(ori_data_b[6])
                        self._nonbonds.append((numa, numb))
                        self._nonbonddatas[(numa, numb)] = ((ra+rb)/0.17818, np.sqrt(epa*epb))
                        break
                    else:
                        raise DataTypeError("can't find van Der Waals data for this atom: {0}".format(str(numb)))
                break
            else:
                raise DataTypeError("can't find van Der Waals data for this atom: {0}".format(str(numa)))
    
    def __gettypes(bondingmap):
        for num, at in enumerate(bondingmap.sites):
            if at[0] == "C":
                Atoms(num, bondingmap)
        for num, at in enumerate(bondingmap.sites):
            if at[0] == "N":
                Atoms(num, bondingmap)
        for num, at in enumerate(bondingmap.sites):
            if at[0] == "O":
                Atoms(num, bondingmap)
        for num, at in enumerate(bondingmap.sites):
            if at[0] == "H":
                Atoms(num, bondingmap)
        for num, at in enumerate(bondingmap.sites):
            if at[0] not in ["C", "N", "O", "H"]:
                Atoms(num, bondingmap)
        return Atoms.gettypelist()[:]


if __name__ == "__main__":
    mol = Molecule()