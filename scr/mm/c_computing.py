import numpy as np
from ctypes import *

lib = None
path = 'C:\\Users\\ARPU\\Documents\\GitHub\\uusrdazoye\\scr\\mm\\C\\mm_c\\Debug\\'

def init():
    global lib
    lib = WinDLL(path + 'mm_c.dll')

def fdihedral_c(num2, array1, array2, array3, array4, dihedraldata):
    ans = np.asarray((0, 0, 0), dtype = 'double')
    p_ans = ans.ctypes.data_as(c_char_p)
    array1, array2, array3, array4, dd = map(lambda x: np.asarray(x).ctypes.data_as(c_char_p), (array1, array2, array3, array4, dihedraldata[1]))
    lib.fdihedral_c(array1, array2, array3, array4, dd, len(dihedraldata[1]), p_ans)
    return ans

def fdihedral_s(num1, array1, array2, array3, array4, dihedraldata):
    ans = np.asarray((0, 0, 0), dtype = 'double')
    p_ans = ans.ctypes.data_as(c_char_p)
    array1, array2, array3, array4, dd = map(lambda x: np.asarray(x).ctypes.data_as(c_char_p), (array1, array2, array3, array4, dihedraldata[1]))
    lib.fdihedral_s(array1, array2, array3, array4, dd, len(dihedraldata[1]), p_ans)
    return ans