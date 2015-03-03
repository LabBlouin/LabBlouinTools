#!/usr/bin/env python

'''
This script is to slice PCA plot into pieces along PC1 coordinate in order to provide principle componet a better
meaning when illustrating clustering display. Also, average will be calculated for each chuck of the plot.
A close model will be selected to represent the cooresponding chuck by comparing to the average using RMSD.
'''

# Date:   Feb 25 2015
# Author: Simiao Lu
# E-mail: simiao.lu@dal.ca

# Import packages #
import sys
from labblouin.PDBnet import PDBstructure as P
from labblouin.Ordination import ORD as O
import numpy as np
from sklearn.decomposition import PCA
from optparse import OptionParser as parser
# End of packages importing #

# Define functions #

def PrepData4PDBnet(prefix):
	traj = P(prefix+'.pdb')
	names = traj.GetModelNames()
	return traj, names

def LoadGMData(prefix):
	coords,labels = [],[]
	with open(prefix+'.gm','r') as F:
		for line in F:
			line=line.strip().strip(';').split(';')
			labels.append(line[0].strip('>'))
			coords.append([float(x) for x in line[1:]])			
	data = np.matrix(coords)
	return data,labels

def PCAcal(prefix,ncomp):
	pc1,pc2,pc3 = [],[],[]
	data,labels = LoadGMData(prefix)
	clf = PCA(ncomp)
	pca = clf.fit_transform(data)
	pcalis = np.array(pca).tolist()
	outf = open('pca_comp3_coords.txt','wb')
	for com in pcalis:
		outf.write('%s %s %s\n'%(com[0],com[1],com[2]))
	outf.close()	
	for m in pcalis:
				pc1.append(m[0])
				pc2.append(m[1])
				pc3.append(m[2])	
	return pc1,pc2,pc3		
	
def PCgenerator(prefix,pc,iv):
	d = {}
	keys = []
	pc1,pc2,pc3 = PCAcal(prefix,ncomp)
	pc_index = int(pc[-1])-1
	pclis = [pc1,pc2,pc3][pc_index]
	for num in range(int(min(pclis))-1,int(max(pclis))+1,int(iv)):
		keys.append(num)
	keys.append(keys[-1]+10)	
	#labels = [[]] * len(keys)
	labels={}
	for k in range(len(keys)):
		d[k] = []
		labels[k]=[]
		for x in range(len(pclis)):
			if keys[k]<pclis[x]<keys[k+1]:
				d[k].append(pclis[x])
				labels[k].append(x)
	return keys,d,labels	

def BelongingClusters(prefix,pc,iv):
	keys,d,labels=PCgenerator(prefix,pc,iv)
	with open(prefix+'.community','r') as F:
		for line in F:
			comp = line.strip().split()
	cluster = {}		
	for k in range(len(keys)):
		cluster[k]=[]
		for i in labels[k]:
			cluster[k].append(comp[i])
	return cluster		
					
def GroupToWrite(prefix,pc,iv):
	keys,d,labels=PCgenerator(prefix,pc,iv)
	s = P()
	traj,names = PrepData4PDBnet(prefix)
	for k in range(len(keys)):
		s = P()
		for i in labels[k]:
			model = traj.GetModel(i)
			s.AddModel(i,model)
		s.WriteFile('%s_%d.pdb'%(pc,k))

# End of functions defining #

# Application #
ncomp = 3
prog = sys.argv[0]
if len(sys.argv) != 4:
	print "Usage: python %s <GMfile_prefix> <interval> <WhichPC>"%(prog)
else:	
	prefix = sys.argv[1]
	iv = sys.argv[2]
	pc = sys.argv[3]	
	GroupToWrite(prefix,pc,iv)