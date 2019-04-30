PT = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 
    'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al':13, 'Si': 14, 'P':15, 'S': 16, 
    'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 
}

class File_IO():

    '''文件IO'''

    def __init__(self, root, handle):
        self.__M = handle
        self.root = root

    def Gaussian_file(self, path):
        inp = open(path, 'r')

        i = 0
        inp_f = {0:[]}
        for line in inp:
            l = line.split()
            if l == []:
                i += 1
                inp_f[i] = []
                continue
            inp_f[i].append(l)
        
        atoms = self.__M.modify_atoms()
        for l in inp_f[2][1:]:
            atoms.append(l[0:1] + [float(ll) for ll in l[1:]])

        bonds = self.__M.modify_bonds()
        for l in inp_f[3]:
            link = []
            for j in range((len(l)-1)//2):
                link.append([ int(l[2*j+1])-1, float(l[2*j+2]) ])
            bonds.append(link)
