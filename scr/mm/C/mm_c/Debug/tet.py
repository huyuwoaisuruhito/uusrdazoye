from ctypes import *
import numpy as np
import os, sys

is_64bits = sys.maxsize > 2**32
if is_64bits:
    print('not suport')
    exit()

CUR_PATH=os.path.dirname(__file__)
dllPath=os.path.join(CUR_PATH,"mm_c.dll")
print (dllPath)
#mydll=ctypes.cdll.LoadLibrary(dllPath)
#print mydll
lib = WinDLL(dllPath)


arr = np.zeros((3,5))
tmp = np.asarray(arr)
rows, cols = tmp.shape
dataptr = tmp.ctypes.data_as(c_char_p)
lib.show_matrix(dataptr, rows, cols)