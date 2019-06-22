from ctypes import *

class ForcesPointer(Structure):
    _fields_ = [('x', c_double), ('y', c_double), ('z', c_double)]