#!/usr/bin/python
'''
this will plot a gm of a 2D shape
'''
import sys
import matplotlib.pyplot as plt
from numpy import mean

if __name__ == '__main__':

	prefix = sys.argv[1] 

	X=[]
	Y=[]
	samples=0
	with open(prefix+'.gm') as F:
		for line in F:
			if line != '':
				X.append([])
				Y.append([])
				bl=line.strip().split(';')
				if any([x == '' for x in bl[1:]]):
					data = bl[1:-1]
				else:
					data = bl[1:]
				for x in range(0,len(data),2):
					X[samples].append(data[x])
				for y in range(1,len(data),2):
					Y[samples].append(data[y])
				samples+=1
	npoints = len(data)/2
	avx=[]
	avy=[]
	for p in range(npoints):
		avx.append([])
		avy.append([])
		temp=[]
		tempy=[]
		for s in range(len(X)):
			temp.append(X[s][p])
			tempy.append(Y[s][p])
		avx[p].append(mean([float(x) for x in temp]))
		avy[p].append(mean([float(y) for y in tempy]))
	
		
		
		
	fig = plt.figure()           
	ax = fig.add_subplot(111)
	#ax.spines['top'].set_color('none')
	#ax.spines['bottom'].set_color('none')
	#ax.spines['left'].set_color('none')
	#ax.spines['right'].set_color('none')
	#ax.xaxis.tick_bottom()
	#ax.xaxis.tick_top()
	#ax.yaxis.tick_left()
	ax.plot(X,Y,'.')
	plt.axis('off')

	for po in range(0,len(avx)):
			ax.annotate(str(po+1),(avx[po][0],avy[po][0]), xycoords='data', xytext=(0, 15), 
				    textcoords='offset points', size=18, 
				    #arrowprops=dict(arrowstyle="fancy",fc="0.6", ec="none", 
				     #               connectionstyle="angle3,angleA=0,angleB=-90")
				     )		            
	plt.savefig('%sPlot.png'%(prefix), dpi=400, format='png')
