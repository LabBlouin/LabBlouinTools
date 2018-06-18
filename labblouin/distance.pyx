__author__ = 'kestrel'

cimport cython
from libc.math cimport sqrt

def dist(A, B):

    cdef double x1 = A[0]
    cdef double y1 = A[1]
    cdef double z1 = A[2]

    cdef double x2 = B[0]
    cdef double y2 = B[1]
    cdef double z2 = B[2]

    cdef double xdif = x1-x2
    cdef double ydif = y1-y2
    cdef double zdif = z1-z2

    cdef int p = 2

    return sqrt(pow(xdif,p) + pow(ydif,p) + pow(zdif,p))
