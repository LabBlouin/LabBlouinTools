import numpy as np

class Candelon_Metiu_Test:
	def __init__(self,observed,nboot=10000):
		'''
		A Distribution-Free Test for Outliers modified from Candelon & Metiu 
		Bundesbank Discussion Paper No /2013
		'''	
		from statsmodels.nonparametric.kde import KDEUnivariate
		from scipy.stats import norm
		self.norm=norm
		self.KDEUnivariate=KDEUnivariate
		self.data=observed
		self.n=len(self.data)
		self.nboot=nboot
		self.Ks={}
		self.bootstrap()
		
	def bootstrap(self):
		size=(self.nboot,len(self.data))
		self.boot = np.random.choice(self.data,size=size,replace=True)
	
	def optimize_lambda(self,data,value,alpha):
		nominal = None
		lalpha=-1
		dist = self.norm(np.mean(data), np.std(data))
		while nominal != alpha:
			v=value*lalpha
			nominal = dist.pdf(v)
			print v,nominal
			if nominal !=0.05:
				lalpha+=0.001
			if lalpha > 2:
				print 'lalpha not found! setting it to 1'
				lalpha=1.0
				return lalpha
			
	def prob(self,data,value,nominal=0.05):
		'''Compute the probability of the values in data being less or equal than value'''
		# set lalpha to the value reported by Hall and York (2001)
		lalpha = 1.1294 #self.optimize_lambda(data, value, nominal)
		c=float(len(np.where(data <= lalpha*value)[0]))
		prob = c/float(len(data))
		return prob >= 1 - nominal
		
	def silverman(self,distribution,hcrit,repeats=1000):
		'''approximation of Silverman's algorithm '''
		Hs=[]
		for i in xrange(repeats):
			# obtain boostrap from optimized distribution
			Mstar=np.random.choice(distribution.support,size=self.nboot)
			# obtain optimized kernel density
			Fstar=self.KDEUnivariate(Mstar)
			Fstar.fit(bw='silverman')
			hstar=Fstar.bw
			Hs.append(hstar)
		return Hs
		
	def BootlierTest(self,k):
		sort = np.array(self.boot)
		sort.sort()
		# remove highest and lowest k
		if k == 0:
			newd = sort
		else:
			newd = sort[:,k:-k]
		# compute the k-trimmed mean of the bth bootstrap
		avYk = newd.mean(axis=1)
		# Compute the difference between the arithmetic mean and the k-trimmed mean
		MTM  = np.array([np.array(sort[b]).mean()-avYk[b] for b in xrange(sort.shape[0])])
		# Obtain the kernel and optimized bandwidth
		KDE=self.KDEUnivariate(MTM)
		KDE.fit(bw='silverman')
		hcrit = KDE.bw
		hstar = self.silverman(KDE, hcrit)
		boolean = self.prob(hstar, hcrit)
		self.Ks[k]=boolean
	
	def outliersID(self):
		'''find the largest subsample which exhibits unimodality'''
		track=[]
		for k in xrange(self.n/2):
			self.BootlierTest(k)
			track.append(self.Ks[k])
			if self.Ks[k] == False:
				return k