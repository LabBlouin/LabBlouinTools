#!/bin/python

''' 
Ordination is a class designed to compute and plot ordination methods such as PCA and MDS.
It is intended as a helper function to PDBnet, but have the functionality to work with 
gm files.

This file is based on Nelle Varoquaux <nelle.varoquaux@gmail.com> code plotmds.py, available
at http://scikit-learn.org/stable/auto_examples/manifold/plot_mds.html, and recomendations in
stackoverflow by Jaime Fernandez (http://numericalrecipes.wordpress.com/)

Dependencies: SKlearn, PDBnet, Scipy, matplotlib

Author: Jose Sergio Hleap
email: jshleap@dal.ca
'''

# importing bit###################################################################################
try:
	import numpy as np
	import scipy.spatial.distance as sp
except:
	print "Dependency Scipy not installed. Please use sudo easy_install scipy, or make sure is installed"

try:
	from sklearn import manifold
	from sklearn.metrics import euclidean_distances
	from sklearn.decomposition import PCA
	from sklearn.lda import LDA
	from sklearn.qda import QDA
except:
	print "Dependency sklearn not installed. Please use sudo easy_install scikit-learn, or make sure is installed"

try:
	from labblouin.PDBnet import PDBstructure as P
except:
	print "Dependency PDBnet not installed or labblouin folder not in pythonpath. Please make sure is set up correctly"

try:
	from mpl_toolkits.mplot3d import Axes3D
	from matplotlib.patches import Ellipse
	import matplotlib.pyplot as plt
	import matplotlib.cm as cm
	from pylab import ion,ioff,draw,plot
except:
	print "Dependency Matplotlib not installed. Please use sudo easy_install matplotlib, or make sure is installed"

from random import shuffle
import pickle
from os.path import isfile
from glob import glob as G
from time import sleep

# END importing bit###############################################################################

# Constants ######################################################################################
colors=['b','r','k','g','y','m','c']
hexa  =['#ED9121', '#EE8262', '#EE1289', '#556B2F', '#FF8C00', '#8B7B8B', '#0000EE', '#EED5D2', 
    '#BA55D3', '#912CEE', '#2F4F4F', '#D15FEE', '#008B8B', '#B23AEE', '#8B7765', '#54FF9F',
    '#8B8386', '#FF4040', '#EEA9B8', '#388E8E', '#6E8B3D', '#33A1C9', '#EE3A8C', '#FF00FF',
    '#436EEE', '#8B864E', '#808000', '#1874CD', '#BCD2EE', '#A9A9A9', '#F4A460', '#FF3030',
    '#FFEBCD', '#B0C4DE', '#00CDCD', '#C0FF3E', '#FFD700', '#8B4513', '#4EEE94', '#CD3278',
    '#00E5EE', '#E3A869', '#CD853F', '#ADD8E6', '#CD2990', '#EEE5DE', '#66CD00', '#7B68EE',
    '#FFA54F', '#A2B5CD', '#BC8F8F', '#8B2323', '#EE30A7', '#EEEED1', '#AEEEEE', '#5E2612',
    '#FF7F00', '#FFC0CB', '#EE3B3B', '#9370DB', '#848484', '#292421', '#CDBA96', '#B4EEB4',
    '#40E0D0', '#8B795E', '#3D9140', '#CDB7B5', '#CAE1FF', '#F0FFFF', '#2E8B57', '#FF6103',
    '#87CEEB', '#CD00CD', '#CDAA7D', '#836FFF', '#EEB4B4', '#8B7355', '#F0E68C', '#CDCDB4',
    '#B4CDCD', '#F0FFF0', '#00EEEE', '#708090', '#9AFF9A', '#FFA07A', '#FFB5C5', '#00688B',
    '#8A3324', '#191970', '#308014', '#FF83FA', '#838B8B', '#808A87', '#00FF7F', '#FFA500',
    '#EEAD0E', '#CD3333', '#4876FF', '#7CCD7C', '#EE5C42', '#AAAAAA', '#DAA520', '#8B3A3A',
    '#FFFAF0', '#B2DFEE', '#00EE76', '#FFFAFA', '#800080', '#C5C1AA', '#EEE685', '#FF3E96',
    '#EE0000', '#FDF5E6', '#EECFA1', '#8DB6CD', '#FF7256', '#7CFC00', '#838B83', '#BF3EFF',
    '#8B6914', '#00CD66', '#A4D3EE', '#00868B', '#8DEEEE', '#8B1C62', '#CDBE70', '#9F79EE', 
    '#C1CDC1', '#CD69C9', '#E0EEEE', '#8B7E66', '#8A2BE2', '#CDCD00', '#97FFFF', '#EEAEEE', 
    '#DC143C', '#CD919E', '#528B8B', '#CD6889', '#E6E6FA', '#E3CF57', '#4B0082', '#FF9912',
    '#F0F8FF', '#FF7F50', '#6CA6CD', '#8B8B83', '#F4F4F4', '#548B54', '#48D1CC', '#C1CDCD', 
    '#E0EEE0', '#3D59AB', '#FFB90F', '#FFD39B', '#8B5A2B', '#9C661F', '#EEE9BF', '#BCEE68',
    '#8EE5EE', '#8B0A50', '#FFF68F', '#EEA2AD', '#CD5B45', '#7FFF00', '#8B8378', '#9BCD9B',
    '#EEE8AA', '#8E8E38', '#668B8B', '#B3EE3A', '#00C78C', '#FFC125', '#8B475D', '#D8BFD8',
    '#FFE4C4', '#96CDCD', '#CDB5CD', '#00C5CD', '#00CED1', '#008B00', '#B8860B', '#1C86EE',
    '#EEC591', '#E066FF', '#B7B7B7', '#DEB887', '#FF6EB4', '#6959CD', '#90EE90', '#8B4789',
    '#EE7AE9', '#8968CD', '#D2B48C', '#FFFFE0', '#CDC9C9', '#BDB76B', '#00C957', '#EEDC82',
    '#3CB371', '#F5FFFA', '#B9D3EE', '#F5F5DC', '#0000CD', '#FF8247', '#EED5B7', '#FFEC8B',
    '#EE7600', '#7171C6', '#8B636C', '#8B814C', '#FFE4B5', '#1E1E1E', '#4F94CD', '#CDAD00', 
    '#CD5555', '#71C671', '#8B7500', '#473C8B', '#B0E0E6', '#FFFF00', '#8B4C39', '#006400',
    '#53868B', '#8B2252', '#FFB6C1', '#63B8FF', '#FFAEB9', '#EE6A50', '#87CEFF', '#87CEFA',
    '#5B5B5B', '#ADFF2F', '#008B45', '#EE4000', '#8A360F', '#8B6969', '#00008B', '#DB7093',
    '#7EC0EE', '#EE799F', '#CD6090', '#C76114', '#8B8682', '#458B74', '#FFF5EE', '#76EE00',
    '#000080', '#228B22', '#8B8B00', '#CD950C', '#EE82EE', '#282828', '#F5DEB3', '#3A5FCD',
    '#00FA9A', '#C67171', '#D1EEEE', '#8B5742', '#8B3E2F', '#CD3700', '#9AC0CD', '#555555',
    '#8B8989', '#EED8AE', '#551A8B', '#778899', '#FFFACD', '#458B00', '#008000', '#FFFFF0',
    '#EEB422', '#5CACEE', '#CD4F39', '#CDC0B0', '#FF7D40', '#8E388E', '#6E7B8B', '#CDC673',
    '#7A378B', '#E0FFFF', '#FFFFFF', '#6C7B8B', '#FFC1C1', '#8B4726', '#515151', '#CD9B1D',
    '#FF6347', '#FF34B3', '#FF0000', '#B0E2FF', '#8B3A62', '#CD5C5C', '#A2CD5A', '#00EE00',
    '#FF6A6A', '#CD6600', '#FFEFDB', '#E9967A', '#EEE9E9', '#7A67EE', '#CD8162', '#00F5FF',
    '#FFEFD5', '#CDAF95', '#00BFFF', '#CDB79E', '#1E90FF', '#EE2C2C', '#8B6508', '#FF7F24',
    '#8FBC8F', '#66CDAA', '#6495ED', '#EEE0E5', '#C1C1C1', '#B22222', '#EE00EE', '#FF82AB',
    '#AB82FF', '#79CDCD', '#7D26CD', '#03A89E', '#8B008B', '#5D478B', '#8B3626', '#808069',
    '#FFE4E1', '#EEDFCC', '#9400D3', '#BFEFFF', '#8B7D6B', '#FF8C69', '#C6E2FF', '#FF4500',
    '#FFE7BA', '#872657', '#808080', '#EE9572', '#CD8500', '#8B5A00', '#9932CC', '#EECBAD',
    '#CD8C95', '#CD1076', '#7D9EC0', '#104E8B', '#8B668B', '#698B22', '#EEE8CD', '#DDA0DD',
    '#4169E1', '#DA70D6', '#DCDCDC', '#68228B', '#CDC8B1', '#000000', '#6B8E23', '#FF69B4',
    '#800000', '#5F9EA0', '#8B4500', '#FCE6C9', '#D3D3D3', '#CDB38B', '#607B8B', '#F08080',
    '#CD9B9B', '#76EEC6', '#FAEBD7', '#68838B', '#EAEAEA', '#7FFFD4', '#C0C0C0', '#EE9A49',
    '#4A708B', '#008080', '#7AC5CD', '#98F5FF', '#8B2500', '#FFF0F5', '#8B8970', '#8B8878',
    '#6A5ACD', '#4682B4', '#EEEEE0', '#27408B', '#00FF00', '#FFDEAD', '#CD2626', '#CD96CD',
    '#9B30FF', '#36648B', '#F8F8FF', '#EEC900', '#EEEE00', '#FFE1FF', '#C1FFC1', '#CDC5BF',
    '#A0522D', '#8B5F65', '#CDC1C5', '#EE7621', '#FFBBFF', '#CD6839', '#698B69', '#BDFCC9',
    '#CD661D', '#FAFAD2', '#CDCDC1', '#FFF8DC', '#B452CD', '#8E8E8E', '#8470FF', '#483D8B',
    '#BBFFFF', '#0000FF', '#EE6AA7', '#EE7942', '#00CD00', '#9ACD32', '#C71585', '#EE9A00',
    '#CAFF70', '#32CD32', '#8B0000', '#B0171F', '#98FB98', '#8B1A1A', '#00B2EE', '#20B2AA',
    '#009ACD', '#A52A2A', '#EE6363', '#FAF0E6', '#8B7D7B', '#9A32CD', '#FF8000', '#7A8B8B',
    '#CD7054', '#9FB6CD', '#CDC9A5', '#D02090', '#00FFFF', '#CD0000', '#43CD80', '#FA8072',
    '#FFDAB9', '#D2691E']
shuffle(hexa)
colors.extend(hexa)
# END Constants ##################################################################################

# Functions bit###################################################################################

class ORD:
	'''
	A class for the most popular ordination methods using PDBnet instaces or gm files.
	'''
	def __init__(self,prefix,data,fastafile=None,n_comp=2):
		self.prefix = prefix
		self.fasta  = fastafile
		self.ncomp  = n_comp

		# vessel for ellipses if created
		self.ellipses={}
		#check if data is PDBnet instance or a filename
		if isinstance(data, P): self.PrepData4PDBnet(data)
		elif isinstance(data,str): self.LoadDataGMfile(data)

		# Center the data
		self.data -= self.data.mean(axis=0)

	def PrepData4PDBnet(self,data):
		'''Load the data to the class assuming is a PDBnet instance file'''
		if data.GetModelNames(): ismodel = True
		else: ismodel = False
		self.labels,self.data = data.gm(self.fasta,ismodel=ismodel,typeof='matrix')


	def LoadDataGMfile(self,data):
		'''Load the data to the class assuming is a GM file'''
		coords , labels = [],[]
		with open(data) as F:
			for line in F:
				bl=line.strip().strip(';').split(';')
				labels.append(bl[0].strip('>'))
				coords.append([float(x) for x in bl[1:]])
		self.data   = np.matrix(coords)
		self.labels = labels

	def dict2array2matrix(self,dict):
		'''
		Giving an upper-triangle distance matrix in a dictionary, returns a distance-like array
		'''
		ar=[]
		keys = sorted(dict.keys())
		for k in sorted(keys):
			for ke in sorted(keys):
				if ke > k:
					ar.append[dict[k][ke]]
		dist = sp.squareform(np.array(ar))
		return dist

	def Store(self):
		self.fit.tofile(self.prefix+'_vecs.temp',sep='\t')
		fout=open(self.prefix+'_infovecs.temp','w')
		fout.write('Number of components:\t%s\nType of ordination:\t%s'%(self.ncomp,self.type))
		fout.close()
		pickle.dump(self, open('%s.current%s.pckl'%(self.prefix,self.type),'wb'))

	def MDS(self,options,typeof='classic',dist=False,groups=None):
		'''
		Perform Multidimensional Scaling wither classic (PCoA) or non-metric.
		If you have the upper triangle of a distance matrix as a dictionary,
		pass the dictionary as dist.
		'''
		# Rotation instance
		self.clf = PCA(n_components=self.ncomp)		

		seed = np.random.RandomState(seed=3)

		if typeof == 'classic': 
			metric = True
			self.type = 'cMDS'
		else: 
			metric = False
			self.type = "nMDS"

		if dist:
			similarities=self.dict2array2matrix(dist)
		else:
			similarities = euclidean_distances(self.data)

		# Initiate multidimensional scaling
		mds = manifold.MDS(n_components=self.ncomp, metric = metric, max_iter=3000, eps=1e-9, 
		                   random_state=seed, dissimilarity="precomputed", n_jobs=-1)

		#fit the data the MDS
		pos = mds.fit(similarities).embedding_
		if typeof != 'classic': pos = mds.fit_transform(similarities, init=pos)

		# Rescale the data
		pos *= np.sqrt((np.array(self.data)** 2).sum()) / np.sqrt((pos ** 2).sum())

		# Rotate the data
		self.fit = self.clf.fit_transform(pos)
		self.Plot(options, groups=groups)
		self.Store()


	def PCA(self,options,groups=None):
		'''Performs a principal coordinate analysis of the data'''
		# Rotation instance
		self.clf = PCA(n_components=self.ncomp)		
		pca = self.clf.fit_transform(self.data)
		self.type='PCA'
		self.fit = pca
		self.explained_variance_ratio = self.clf.explained_variance_ratio_
		self.Plot(options,groups=groups)
		self.Store()
	
	def loopoverpcs(self,prefix,comp,options):
		samemarker = ['b.','r.','k.','g.','y.','m.','c.']
		for z in zip([0,1,0],[1,2,2]):
			fig2 = plt.figure(dpi=options.dpi)
			ax2  = fig2.add_subplot(111)
			ax2.spines['top'].set_color('none')
			ax2.xaxis.tick_bottom()
			ax2.spines['right'].set_color('none')
			ax2.yaxis.tick_left()
			for i in xrange(len(comp)):
				point = comp[i,:]
				color = str(i/float(len(comp)))				
				ax2.plot(point[z[0]],point[z[1]],c=color,marker=samemarker,markersize=options.markersize,linestyle='None',alpha=a)
				ax.set_xlabel('PC %d (%.2f%%)'%(z[0],self.explained_variance_ratio[0]*100), fontsize=fontsize)
				ax.set_ylabel('PC %d (%.2f%%)'%(z[1],self.explained_variance_ratio[1]*100), fontsize=fontsize)				
			fig2.savefig(prefix+'PC%d.PC%d.%s'%(z[0],z[1],options.of))
			fig2 = plt.figure(dpi=options.dpi)	
				

	def Plot(self,options,groups=None):
		'''
		Plot the components from an ordination method of the class ORD. If the number of components
		is greater than 3, it will plot the first three components. Components has to be a n x k 
		numpy array of eigenvectors, where n is the observations/individuals and k the components.
		The option groups allow to pass a list (of the same lenght of the arrar, that is a lenght of n).
		'''
		dpi = options.dpi
		fontsize = options.textsize

		components=self.fit
		samemarker = ['b.','r.','k.','g.','y.','m.','c.']
		markers = ['k.','b+','g*','r.','c+','m*','y.','k+','b*','g.','r+','c*','m.','y+','k*','b.',
		           'g+','r*','c.','m+','y*']	
		if options.interactive:
			ioff()
		if components.shape[1] >= 3: 
			dim = 3
			fig = plt.figure(dpi=dpi)
			#ax  = fig.gca(projection='3d')
			ax = Axes3D(fig)
			

		else: 
			dim = 2
			fig = plt.figure(dpi=dpi)
			ax  = fig.add_subplot(111)
			ax.spines['top'].set_color('none')
			ax.xaxis.tick_bottom()
			ax.spines['right'].set_color('none')
			ax.yaxis.tick_left()

		comp = components[:, 0:dim]

		#if MD:
		#	colormap=cm.gray
		if options.samemarker:
			markers = samemarker

		if options.interactive:
			draw()
		if groups != None and dim == 2:			
			if len(set(groups)) <= len(markers):
				g=sorted(list(set(groups)))
				for gr in g:
					c  = markers[g.index(gr)][0]
					m  = markers[g.index(gr)][1]
					ms = options.markersize + 10
					#d[gr]=markers[list(g).index(gr)]	
					x = comp[np.where(groups == gr),0]
					y = comp[np.where(groups == gr),1]
					ax.scatter(x, y, c= c, marker=m,label=gr,s=ms)
					if options.numbered:
						gg = len(comp[np.where(groups == gr)])	
						for num in range(gg):
							ax.annotate('%d'%(num),xy=(list(x[0])[num],list(y[0])[num]),xycoords='data')
					#markersize=options.markersize,linestyle='None',label=gr)#c=color,marker=marker
						

		elif options.interactive:
			ion()
			count=-1
			step = range(5,len(comp),50)
			for i in xrange(len(comp)):
				count += 1
				point = comp[i,:]
				if options.MD and i == 0:
					start=point
					color='red'
					marker=r'$\circlearrowleft$'
					a = None
				elif options.MD and i == len(comp)-1:
					last=point
					color='blue'
					marker=r'$\lambda$'
					a = None
				else:
					marker='o'
					color = str(i/float(len(comp)))				
					a = 0.5

				ax.plot(point[0],point[1],c=color,marker=marker,markersize=options.markersize,linestyle='None',alpha=a)
				if options.interactive:
					plot(point[0],point[1],c=color,marker=marker,markersize=options.markersize,linestyle='None',alpha=a)
					if count < 10:	
						sleep(1)
						draw()
					elif i in step:
						draw()				



			if options.MD:
				ax.plot(start[0],start[1],c="red",marker=r'$\circlearrowleft$',markersize=options.markersize,linestyle='None')
				ax.plot(last[0] , last[1],c="blue",marker=r'$\lambda$',markersize=options.markersize,linestyle='None')

			#ax.scatter(comp[:,0],comp[:,1],c=color,cmap = colormap, s=15)
			
			groups=[None]*len(self.labels)
		elif groups != None and dim == 3:
			g=sorted(list(set(groups)))
			for gr in g:
				c  = markers[g.index(gr)][0]
				m  = markers[g.index(gr)][1]
				ms = options.markersize + 10
				x = comp[np.where(groups == gr),0]
				y = comp[np.where(groups == gr),1]
				z = comp[np.where(groups == gr),2]
				ax.scatter(x,y,z,color=c,marker=m,label=gr,s=ms)
				if options.rotate > 0:
						for ii in xrange(0,360,options.rotate):
							ax.view_init(elev=0, azim=ii)
							ax.set_xlabel('PC 1 (%.2f%%)'%(self.explained_variance_ratio[0]*100), fontsize=fontsize)
							ax.set_ylabel('PC 2 (%.2f%%)'%(self.explained_variance_ratio[1]*100), fontsize=fontsize)
							ax.set_zlabel('PC 3 (%.2f%%)'%(self.explained_variance_ratio[2]*100), fontsize=fontsize)
							fig.savefig(self.prefix + "rotation%d"%(ii) + ".%s"%(options.of),dpi=dpi)				
					
		fout=open('%s.equivalences'%self.prefix,'w')
		if not options.MD:
			for l in range(len(self.labels)):
				ax.annotate(l,comp[l,]+0.1,fontsize=fontsize)
				fout.write('%s\t%s\t%s'%(l,self.labels[l],groups[l]))
					
		if self.type == 'PCA':
			ax.set_xlabel('PC 1 (%.2f%%)'%(self.explained_variance_ratio[0]*100), fontsize=fontsize)
			ax.set_ylabel('PC 2 (%.2f%%)'%(self.explained_variance_ratio[1]*100), fontsize=fontsize)
			if dim >= 3:
				ax.set_zlabel('PC 3 (%.2f%%)'%(self.explained_variance_ratio[2]*100), fontsize=fontsize)
				ax.view_init(30, 45)
		else:
			ax.set_xlabel('Axis 1', fontsize=fontsize)
			ax.set_ylabel('Axis 2', fontsize=fontsize)
			if dim >= 3:
				ax.set_zlabel('Axis 3', fontsize=fontsize)
				ax.view_init(30, 45)
		if groups[0]:
			#handles, labels = ax.get_legend_handles_labels()
			ax.legend(loc=0, fancybox=True, shadow=True, fontsize='small')

		if dim == 2:
			fig.tight_layout()
		#plt.show()
				
		else:
			fig.savefig(self.prefix+'_%s.%s'%(self.type,options.of), dpi=dpi)
		fout.close()		

		if options.MD:
			initx,inity = comp[:,0][0],comp[:,1][0]
			lastx,lasty = comp[:,0][-1],comp[:,1][-1]
			ax.annotate('Initial conformation', xy=(initx, inity), 
			            xytext= (float(initx)+10, float(inity)+20),
			            arrowprops=dict(facecolor='blue', shrink=0.05, frac=0.15))

			ax.annotate('Final conformation', xy=(lastx, lasty),
			            xytext=(float(lastx)+10, float(lasty)+40),
			            arrowprops=dict(facecolor='blue', shrink=0.05, frac=0.15))
			fig.savefig(self.prefix+'_%s_startStop.%s'%(self.type,options.of), dpi=dpi)
		

	def PlotXDA(self,membership,options,group_labels=None,):
		'''
		Plots a Linear Discriminant Analysis (LDA) or a Quadratic Discriminan Analysis (QDA) with
		confidence ellipses at std (standard deviations)
		'''
		#options:
		stds=options.std
		ellipses=options.ellipses
		dpi=options.dpi
		typeof=self.type
		fontsize=options.textsize
		MD=options.MD
		legend = options.legend
		numbered = options.numbered

		if group_labels:
			target_names = group_labels
		else:
			target_names=sorted(list(set(membership)))	
		fig = plt.figure()
		ax = fig.add_subplot(111)
		for c, i, target_name in zip(colors, list(set(membership)), target_names):
			sub = self.fit[membership == i, ]
			x=sub[:,0]
			Y=sub[:,1]
			if ellipses:
				for j in xrange(1, stds+1):
					ax.add_artist(self.ellipses[i][j])
				ax.scatter(x, Y, c=c, label=target_name)
				if numbered:
					for num in range(len(sub)):
						ax.annotate('%d'%(num),xy=(list(x)[num],list(Y)[num]),xycoords='data')				
				if legend:
					plt.legend(target_names,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
		ax.spines['top'].set_color('none')
		ax.xaxis.tick_bottom()
		ax.spines['right'].set_color('none')
		ax.yaxis.tick_left()
	
		if typeof == 'LDA':
			ax.set_xlabel('LD 1', fontsize=fontsize)
			ax.set_ylabel('LD 2', fontsize=fontsize)			
		else:
			ax.set_xlabel('QD 1', fontsize=fontsize)
			ax.set_ylabel('QD 2', fontsize=fontsize)

		if MD:
			initx,inity = self.fit[:,0][0],self.fit[:,1][0]
			lastx,lasty = self.fit[:,0][-1],self.fit[:,1][-1]
			ax.annotate('Initial conformation', xy=(initx, inity), 
			            xytext= (float(initx)+10, float(inity)+20),
			            arrowprops=dict(facecolor='blue', shrink=0.05, frac=0.15))

			ax.annotate('Final conformation', xy=(lastx, lasty),
			            xytext=(float(lastx)+10, float(lasty)+40),
			            arrowprops=dict(facecolor='blue', shrink=0.05, frac=0.15))		
		plt.show()
		fig.savefig(self.prefix+'_%s.%s'%(typeof,options.of), dpi=dpi)

	def ellipse(self, singlegroupx,singlegroupy,std=2,color='k'):
		''' 
		Create an ellipse given points with x coordinates in singlegroupx and singlegroupy
		'''
		cov        = np.cov(singlegroupx, singlegroupy)
		lambda_, v = np.linalg.eig(cov)
		lambda_    = np.sqrt(lambda_)
		width      = lambda_[0]*std*2
		height     = lambda_[1]*std*2
		centerx    = np.mean(singlegroupx)
		centery    = np.mean(singlegroupy)
		angle      = np.rad2deg(np.arccos(v[0, 0]))
		ell        = Ellipse(xy=(centerx,centery), width=width , height=height , angle=angle, 
		                     color=color)
		ell.set_facecolor('none')

		return ell

	def pointsInEllipse(self,Xs,Ys,ellipse):
		'''
		Tests which set of points are within the boundaries defined by the ellipse. The set of points
		are defined in two arrays Xs and Ys for the x-coordinates and y-coordinates respectively
		'''
		inside = []
		outside= []
		for i in range(len(Xs)):
			point = (Xs[i],Ys[i])
			if ellipse.contains_point(point):
				inside.append(point)
			else:
				outside.append(point)
		return inside, outside

	def getEllipses(self,stds):
		'''will populate the ellipses attribute'''

		for mem,col in zip(list(set(membership)),colors):
			self.ellipses[mem]={}
			sub = self.fit[membership == mem, ]
			Xs=sub[:,0]
			Ys=sub[:,1]
			for j in xrange(1, stds+1):
				self.ellipses[mem][j]=self.ellipse(Xs,Ys,std=j,color=col)

	def LDA(self,membership,options,group_labels=None):
		'''
		Perform a Linear discriminant analysis of the data and plot it. Membership must be an array
		of integers of the same lenght of the number of observations in the data. 
		'''
		self.type='LDA'
		lda = LDA(n_components=2)
		self.fit = lda.fit(self.data, membership).transform(self.data)
		if options.ellipses:
			self.getEllipses(options.std)
		self.PlotXDA(membership, options, group_labels=group_labels)
		self.Store()

	def QDA(self,membership, options, group_labels=None):
		self.type = 'QDA'
		qda = QDA()
		self.fit = qda.fit(self.data, membership).predict(self.data)
		if options.ellipses:
			self.getEllipses(options.std)
		self.PlotXDA(membership,options,group_labels=group_labels)
		self.Store()
# END Functions bit###############################################################################
if __name__ == "__main__":
	import optparse,sys
	opts = optparse.OptionParser(usage='%prog <prefix> <GM file> [options]')
	opts.add_option('-s','--std',dest='std', action="store", default=3,type=int, help='Standard '\
	                'deviations to plot ellipses. Has no effect for PCA or any of the MDS.'\
	                'Default:3')
	opts.add_option('-m','--mem',dest='mem', action="store", default=None, help='File name for '\
	                'the membership vector if any of the discriminant analysis is usde. The file'\
	                ' is a space sepatated format, with the groupings')	
	opts.add_option('-t','--type',dest='type', action="store", default='All', help='Type of '\
	                'analysis. Available: classic MDS (PCoA), non-metric MDS (nMDS), Linear'\
	                ' Discriminan Analysis (LDA), Quadratic discriminant analysis (QDA) or PCA.'\
	                ' You can pass All if you want to do all of the above. Default:All')
	opts.add_option('-M','--MD',dest='MD', action="store_true", default=False, help='If too many '\
	                'snapshots in a molecular dynamic simulation are used, this option will stop '\
	                'printing the labels, and will point at the extremes of the trajectory')
	opts.add_option('-d','--dpi',dest='dpi', action="store", default=300, type=int,
	                help='Modify the dpi of the plots. Default: 300')	
	opts.add_option('-c','--ncomp',dest='ncomp', action="store", default=2, help='Modify the number'\
	                ' of components to be estimated. Default: 2',type=int)		
	opts.add_option('-g','--groups',dest='groups', action="store", default=None, help='Provide a '\
	                'grouping vector to color plots by. This file is a space sepatated format, just'\
	                ' as the membership file. Default=No')	
	opts.add_option('-T','--textsize',dest='textsize', action="store", default=10, type=int,
	                help='Modify the text size in the plot. Default=10')
	opts.add_option('-e','--ellipses',dest='ellipses', action="store_false", default=True, help='Turn off '\
	                'the ellipses drawing.')	
	opts.add_option('-i','--interactive',dest='interactive', action="store_true", default=False, 
	                help='Whether or not to do interactive plotting in such a way that the trajectory'\
	                ' can be seen. This can be really slow if too many points are being drawn. Default = False')		
	opts.add_option('--ms','--markersize',dest='markersize', action="store", default=8, type=int,
	                help='the size of the marker in the plot. Default = 8')
	opts.add_option('-o','--outputformat',dest='of', action="store", default='pdf', help='Select the format of'\
	                ' the plot. pgf: LaTeX PGF Figure, svgz: Scalable Vector Graphics , tiff: Tagged Image File Format'\
	                ', jpg: Joint Photographic Experts Group, raw: Raw RGBA bitmap, jpeg: Joint Photographic Experts Group'\
	                'png: Portable Network Graphics, ps: Postscript, svg: Scalable Vector Graphics, eps: Encapsulated Postscript'\
	                ', rgba: Raw RGBA bitmap, pdf: Portable Document Format, tif: Tagged Image File Format. Default = pdf')		
	opts.add_option('-l','--legend',dest='legend',action='store_true',default=False,help='Indicate if you want to add a legend to the '\
	                'XDA plot. Default = False')
	opts.add_option('--sm', '--samemarker', dest='samemarker', action='store_true', default=False,help='This feature is avaible if you want to plot XDA with'\
	                ' the same marker (dot) but with different colors when # of clusters <= 7. Also, the colors match with a LDA plot with the same data provided.'\
	                ' Default = False.')
	opts.add_option('-n','--numbered',dest='numbered',action='store_true',default=False, help='This feature is available if you want'\
	                ' the PCA or XDA plot with model names indicated. Default = False.')
	opts.add_option('-r','--rotate',dest='rotate',action='store',default=0,help='This option helps you'\
	                ' create a series of snapshots of rotating the 3D PCA plot. Please indicate the stepsize (from 0 to 360 deg).'\
	                ' Default = 0',type=int)
	opts.add_option('-L','--loop',dest='loop',action='store_true',default=False,help='This feature is available to'\
	                ' plot PC1 vs PC2, PC2 vs PC3 and PC1 vs PC3. Default=False')

	options, args = opts.parse_args()
	ncomp = int(options.ncomp)
	if options.type == 'All' and not options.mem:
		print 'You choose to use all the functions but did not provide a membership file.'
		sys.exit()
	if options.mem:
		membership = np.array(open(options.mem).read().strip().split())
	if options.groups:
		groups = np.array(open(options.groups).read().strip().split())
	else:
		groups=None
	if not G('*.pckl'):
		O = ORD(args[0], args[1],n_comp=ncomp)
		if options.type == 'All' or options.type == 'PCoA': O.MDS(options,groups=groups)
		if options.type == 'All' or options.type == 'nMDS': O.MDS(options,typeof='non-metric',groups=groups)
		if options.type == 'All' or options.type == 'PCA' : O.PCA(options,groups=groups)
		if options.type == 'All' or options.type == 'LDA' : O.LDA(membership,options,group_labels=groups)
		if options.type == 'All' or options.type == 'QDA' : O.QDA(membership,options,group_labels=groups)
	else:
		for i in G('*.pckl'):
			O = pickle.load(open(i,'rb'))
			if 'DA' in i:	
				O.PlotXDA(membership,options)
			else:
				O.Plot(options,groups=groups)
		