'''
Helper function to do GPA with ABeRMuSA
'''
import rpy2.robjects as robj
from rpy2.robjects.packages import importr
import numpy as np

def array2R(arr):
	'''
	Input a numpy array into R
	
	:param arr: a 2d numpy array with the data.
	:type arr: :class `numpy.matrix`
	:returns: an rp2 matrix object
	'''
	nr, nc = arr.shape
	xvec = robj.FloatVector(arr.reshape((arr.size)))
	xr = robj.r.matrix(xvec, nrow=nr, ncol=nc,byrow=True)
	return xr

def Rmatrix2Rarray(matrix,dim=3):
	'''
	Convert an rpy matrix into a k x m x n array, where k is the number of landmarks, m the number
	of dimensions, and n the number of observations
	
	:param matrix: an rpy2 matrix of dimensions km x n
	:type matrix: :class `rpy2.robjects.vectors.Matrix`
	:returns: a k x m x n rpy2 array
	'''
	## Convert GMfile into 'shapes'-readable array
	## A is the read.table variable of a GMfile
	## s is the sample size
	## k is the number of landmarks
	## d is the number of dimensions (2D,2;3D,3)
	robj.r('arr <- function(A,d){s<-dim(A)[1]; k<-(dim(A)[2]/d);'\
	  'Data <- as.numeric(A);arr<-array(Data,c(s, d, k)); arr<-aperm(arr,c(3,2,1)); arr}')
	
	return robj.r['arr'](matrix,dim)

def shapesGPA(arr,scale=False):
	'''
	Perform the generilized procrustes superimpostion
	
	:param arr: a 2d numpy array with the data.
	:type arr: :class `rpy2.robjects.vectors.Array`
	:param scale: wether to scale the shape or not (default).
	:type scale: boolean
	:returns: the rotated version on the array
	'''
	shapes = importr('shapes')
	gpa = shapes.procGPA(arr,scale=scale,distances=False, pcaoutput=False)
	return gpa[gpa.names.index('rotated')]

def A2GM(array,dim=3):
	'''
	Convert 'shapes'-readable array into GM-type. Return a numpy array
	
	:param array: a a k x m x n array, where k is the number of landmarks, m 
	the number of dimensions, and n the number of observations.
	:type arrray: :class `rpy2.robjects.vectors.Array`
	:param dim: dimensions
	:type dim: integer
	'''
	## Convert 'shapes'-readable array into GMfile 
	## A is an array Shapes-type
	## d is the dimensions
	robj.r(
	    'A2GM<- function(A,d,rownames=NULL){m<-matrix(NA,dim(A)[3],dim(A)[1]'\
	    '*d); for (i in 1:dim(A)[3]){ for (j in 1:d){ m[i,seq(j,dim(m)[2],d)]'\
	    '<-A[,j,i]}}; row.names(m)<-rownames; return(m)}')	
	mat = robj.r['A2GM'](array,dim)
	return np.array(mat)

def plot_shapes(arr):
	'''
	given a n (number of observation) x k (number of landmarks) x m (number of dimensions), plot
	the superimposition or original arrangement
	
	:param arr: Input (n x k x m)  2d numpy array, where k is the number of points, m is the number
	of dimensions, and n is the sample size.
	:type arr: :class `numpy.array`
	'''
	n,k,m = arr.shape
	fig=plt.figure()
	for i in xrange(n):
		# Get every point in tuple form
		A = tuple(arr[i,:,:])
		for j in xrange(k):
			plt.scatter(A[j][0],A[j][1])
	plt.show()
	
if __name__ == "__main__":
	# for debugging
	import sys
	import matplotlib.pylab as plt

	fn = sys.argv[1]
	dim = int(sys.argv[2])
	scale = bool(sys.argv[3])
	matrix=[]
	with open(fn) as F:
		for line in F:
			bl=line.strip().strip(';').split(';')
			matrix.append([float(x) for x in bl[1:]])
	matrix = np.asarray(matrix)
	row, col = matrix.shape
	k = col/dim
	arr = Rmatrix2Rarray(array2R(matrix), dim=dim)
	rot = shapesGPA(arr,scale=scale)
	mat = A2GM(rot, dim=dim)
	plot_shapes(mat.reshape((row,k,dim)))