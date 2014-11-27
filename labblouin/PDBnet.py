#!/bin/python

""" PDBnet is a collection of Python objects intended to model and contain
PDB protein data.

PDBnet Copyright (C) 2012 Christian Blouin
Contributions by Alex Safatli and Jose Sergio Hleap
Major Refactoring (2014) done by Alex Safatli and Jose Sergio Hleap

E-mail: cblouin@cs.dal.ca, safatli@cs.dal.ca, jshleap@dal.ca
Dependencies: Scipy, BioPython, FASTAnet (contained in LabBlouinTools) """

from mmap import mmap as memoryMappedFile
import numpy as np
import sys, FASTAnet, math
from math import sqrt
from os.path import getsize as sizeOfFile
from Bio.SeqUtils.ProtParam import ProteinAnalysis as PA
from scipy.spatial.distance import cdist

# Metadata

__author__  = ['Christian Blouin','Jose Sergio Hleap','Alexander Safatli']
__version__ = '1.0.5'

# Constants

PDB_LARGE_FILE_SIZE = 100000000 #bytes

# Classes

class PDBatom(object):

	""" ATOM in a PDB protein structure. """

	__slots__ = ('name','serial','x','y','z','occupancy','tempFactor',
	             'charge','symbol','parent')

	def __init__(self,serial,name,x,y,z,oc,b,symbol,charge):

		""" Construct an atom in a PDB protein structure. """

		self.name       = name
		self.serial     = serial
		self.x          = x
		self.y          = y
		self.z          = z
		self.occupancy  = oc
		self.tempFactor = b
		self.charge     = charge
		self.symbol     = symbol
		self.parent     = None

	def fixname(self):

		""" Ensures the name of this atom fits within 4 characters. """

		if len(self.name) == 4: return self.name
		else: return ' '+self.name.ljust(3)

	def __str__(self):

		""" Express this atom as a single line in a PDB file (see PDB format). """

		return "ATOM  %5s %4s %3s %1s%4s    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s" \
		       % (self.serial,self.fixname(),self.parent.name,self.parent.chain,
		          self.parent.index,self.x,self.y,self.z,self.occupancy,
		          self.tempFactor,self.symbol,self.charge)

	def GetPosition(self): 

		""" Get the 3-dimensional coordinate of this atom in space. """

		return (self.x,self.y,self.z)

	def DistanceTo(self, atom):

		""" Acquire the distance from this atom to another atom. """

		return ((self.x-atom.x)**2+(self.y-atom.y)**2+(self.z-atom.z)**2)**0.5

class PDBterminator(PDBatom):

	""" A placeholder class that represents a terminating ATOM-like line in the PDB file. """

	__slots__ = ('name','serial','x','y','z','occupancy','tempFactor',
	             'charge','symbol','parent','lastatom','lastresname',
	             'lastreschain','lastresind')

	def __init__(self,chaininst):

		self.lastatom = chaininst.GetAtoms()[-1]
		self.lastresname = self.lastatom.parent.name
		self.lastreschain =  self.lastatom.parent.chain
		self.lastresind = self.lastatom.parent.index
		newserial = str(int(self.lastatom.serial) + 1)
		super(PDBterminator,self).__init__(newserial,'',0,0,0,0,0,'','')

	def __str__(self):

		return 'TER   %5s %4s %3s %1s%4s' % (
		    self.serial,'',self.lastresname,self.lastreschain,self.lastresind)

class PDBresidue:

	""" A residue (collection of ATOM fields) in a PDB protein structure. """

	__slots__ = ('index','name','atoms','atomsOrdered','contacts','centroid',
	             'chain','model')

	def __init__(self,index=None,name=''):

		""" Construct a residue (collection of atoms) in a PDB protein structure. """

		self.index = index
		self.name = name
		self.atoms = {}
		self.atomsOrdered = []
		self.contacts = []
		self.centroid = None
		self.chain = ''
		self.model = None

	def GetAtoms(self): return [x for x in self.atomsOrdered]

	def AddAtom(self, atom):

		""" Add a PDBatom structure to this residue. """

		self.atoms[atom.name.strip()] = atom
		self.atomsOrdered.append(atom)
		atom.parent = self

	def GetCA(self):

		""" Get the alpha-carbon found in this residue as a PDBatom. """

		return self.atoms['CA']

	def Centroid(self):

		""" Calculate the centroid of this residue. Return this as a PDBatom. """

		x = 0.0
		y = 0.0
		z = 0.0

		for atom in self.atoms:
			if atom in ['C','N','O']: continue
			x += self.atoms[atom].x
			y += self.atoms[atom].y
			z += self.atoms[atom].z

		div = float(len(self.atoms)-3)

		x /= div
		y /= div
		z /= div

		self.centroid = PDBatom(None,'centroid',x,y,z,0,0,'','')
		self.centroid.parent = self
		return self.centroid

	def InContactWith(self,other,thres=4.5):

		""" Determine if in contact with another residue. """

		foundmatch = False
		for atomA in self.atoms:
			if atomA in ['C', 'N', 'O']:
				continue
			atomA = self.atoms[atomA]

			for atomB in other.atoms:
				if atomB in ['C', 'N', 'O']:
					continue   
				atomB = other.atoms[atomB]

				if atomA.DistanceTo(atomB) <= thres:
					foundmatch = True

				if foundmatch:
					break
			if foundmatch: 
				break     

		return foundmatch        

	def __int__(self): return int(self.index)
	def __str__(self): return '\n'.join([str(x) for x in self.atomsOrdered])

class PDBchain(object):

	""" A PDB chain (collection of protein residues). """

	__slots__ = ('name','structure','residues','indices','parent')

	def __init__(self,name):

		self.name       = name
		self.residues   = list()
		self.indices    = list()
		self.structure  = None
		self.parent     = None

	def GetResidues(self): 

		""" Acquire all of the residues in this chain. Note that if this
		is a large file, this will force the object to load all of these
		residues into memory from the file. """

		return [self.GetResidueByIndex(i) for i in self.GetIndices()]

	def IterResidues(self):

		""" Iteratively yield all of the residues in this chain. Note that if this
		is a large file, this will force the object to load all of these
		residues into memory from the file. """

		for i in self.GetIndices(): yield self.GetResidueByIndex(i)

	def GetResidueByIndex(self,index): 

		""" Alias for acquiring a residue from this instance using [] operator. """

		return self.__getitem__(index)

	def GetAtoms(self): 

		""" Get all comprising atoms in this chain in order of residues. """

		atoms = []
		for x in self.IterResidues():
			atoms.extend(x.GetAtoms())
		return atoms

	def GetIndices(self): 

		""" Acquire all of the indices for residues present in this chain. """

		return [int(x) for x in self.indices]    

	def AddIndexOfResidue(self,index):

		""" Add an index for a residue to this chain. Will generate own residue
		information if called upon (forces lazy evaluation). """

		self.indices.append(str(index))

	def AddResidue(self,resid):

		""" Add a residue to this chain. """

		resid.chain = self.name
		self.residues.append(resid)
		self.indices.append(str(resid.index))

	def AddResidueByIndex(self,index):

		""" Add a residue to this chain by its index in the PDB. The residue
		object will automatically be constructed from the file. """

		self.indices.append(str(index))
		resid = self.GetResidueByIndex(index)
		self.residues.append(resid)

	def RemoveResidue(self,resid):

		""" Remove a residue from this chain. """

		return self.pop(resid)

	def ContactMap(self,thres=4.5):

		""" Compute the contact map of this chain. """

		done       = []
		contactmap = []
		for rA in self:
			for rB in self:
				resA,resB = self[rA],self[rB]
				if resA == resB or (resA,resB) in done: continue  
				elif resA.InContactWith(resB,thres):
					if not (int(resA),int(resB)) in contactmap:
						contactmap.append((int(resA),int(resB)))
						contactmap.append((int(resB),int(resA)))
				done.append((resA,resB))
				done.append((resB,resA))
		return sorted(contactmap, key=lambda e: (e[0], e[1]))

	def GetPrimaryPropertiesFromBioPython(self):

		""" Use BioPython to populate this class with attributes. """

		props=PA(self.AsFASTA())
		self.amino_acid_composition= props.get_amino_acids_percent()
		self.molecular_weight= props.molecular_weight()
		self.aromaticity= props.aromaticity()
		self.instability_index= props.instability_index()
		self.flexibility= props.flexibility()
		self.isoelectric_point=props.isoelectric_point()
		self.secondary_structure_fraction=props.secondary_structure_fraction()        
		return props

	def AsFASTA(self):

		""" Return the string representing this chain's FASTA sequence. """

		return ''.join([aa[item.name] for item in self.IterResidues()])

	def WriteAsPDB(self,filename):

		""" Write this single chain to a file as a PDB. """

		fh = open(filename,'w')
		fh.write(str(self))
		fh.close()        

	def SortByNumericalIndices(self):

		""" Sort all internal items by a numerical index. """

		if not self.structure.handle.isLargeFile():
			self.residues = sorted(self.residues,key=lambda d: int(d.index))
		self.indices = sorted(self.indices, key=lambda d: int(d))

	def pop(self,s):

		""" Pop a residue. """

		res = self.GetResidueByIndex(s)
		if not self.structure.handle.isLargeFile():
			self.residues.remove(res)
		self.indices.remove(str(s))
		return res

	def update(self,o):

		""" Add residues from another chain. """

		for res in o: self.AddResidue(res)

	def __str__(self):

		get = lambda d: self.GetResidueByIndex(d)
		out = '\n'.join([str(get(x)) for x in self.indices])
		out += '\n%s\n' % (str(PDBterminator(self)))
		return out

	def __getitem__(self,i):

		""" Force the object to load the residue if not present
		in this structure. """

		handle = self.structure.handle
		if str(i) in self.indices:
			if (handle and handle.isLargeFile()) or (len(self.indices) > len(self.residues)):
				if self.parent == None:
					resid = handle.readResidue(self.name,str(i))
				else:
					resid = handle.readResidue(self.name,str(i),self.parent.name)
				resid.chain = self.name
				return resid
			return self.residues[self.indices.index(str(i))]
		else: return None

	def __len__(self): return len(self.indices)

	def __iter__(self):
		for ind in self.indices: yield int(ind)

class PDBmodel(PDBchain):

	""" A PDB model (a special kind of chain). """

	__slots__ = ('name','structure','residues','indices','chains','chainNames')

	def __init__(self,name):

		if type(name) != int: raise TypeError('Model names should be integers.')
		elif name < 0: raise ValueError('Model name should be zero or a positive integer.')
		super(PDBmodel,self).__init__(name)
		self.chains     = list()
		self.chainNames = list()

	def GetResidues(self): 

		""" Acquire all of the residues in this model. Note that if this
		is a large file, this will force the object to load all of these
		residues into memory from the file. """

		this = [self.GetResidueByIndex(i) for i in self.GetIndices()]
		for it in self.chains: this += it.GetResidues()
		return this

	def IterResidues(self):

		""" Iteratively yield all of the residues in this model. Note that if this
		is a large file, this will force the object to load all of these
		residues into memory from the file. """

		for i in self.GetIndices(): yield self.GetResidueByIndex(i)
		for i in self.chains:
			for j in i.IterResidues(): yield j

	def GetChains(self): return self.chains

	def GetChain(self,name):
		return self.GetChainByName(name)

	def GetChainByName(self,name):
		for x in xrange(len(self.chainNames)):
			if self.chainNames[x] == name: return self.chains[x]

	def GetChainNames(self): return self.chainNames

	def NewChain(self,ch):

		""" Create a new PDBchain, add it to this model, and return it. """

		c = PDBchain(ch)
		self.AddChain(c)
		return c

	def AddChain(self,chain):

		""" Add a PDBchain (chain) to this model. """

		chain.structure = self.structure
		chain.parent    = self
		self.chains.append(chain)
		self.chainNames.append(chain.name)

	def AddResidue(self,resid):

		""" Add a residue to this model. """

		resid.model = self.name
		self.residues.append(resid)
		self.indices.append(str(resid.index))

	def __getitem__(self, i):

		handle = self.structure.handle
		if str(i) in self.indices:
			if (handle and handle.isLargeFile()) or len(self.indices) > len(self.residues):
				resid = handle.readResidue('',str(i),self.name)
				resid.model = self.name
				return resid
			return self.residues[self.indices.index(str(i))]
		else: return None		

	def __str__(self):

		model = 'MODEL%9s\n' % (str(self.name)) + super(
		    PDBmodel,self).__str__()
		for ch in self.chains: model += str(ch)
		return model + 'ENDMDL\n'

class PDBstructure(object):

	""" A PDB protein structure (a collection of chains/models). """

	__slots__ = ('chains','orderofchains','models','orderofmodels','remarks',
	             'filepath','organism','taxid','mutation','contactmap','handle',
	             'ismodel')

	def __init__(self, filein=''):

		# Attributes.
		self.filepath      = filein
		self.handle        = None
		self.organism      = None
		self.taxid         = None
		self.mutation      = False
		self.contactmap    = list()
		self.orderofchains = list()
		self.orderofmodels = list()
		self.remarks       = list()
		self.chains        = dict()
		self.models        = dict()

		# Read a file?
		if filein != '':
			self.ReadFile(filein)

	# Accesssors/Mutators.

	def GetChainNames(self): return self.orderofchains
	def GetModelNames(self): return self.orderofmodels

	def GetChain(self,ch): 

		""" Get a chain by name. """

		if ch in self.chains: return self.chains[ch]
		return None

	def GetModel(self,mod):

		""" Get a model by name. """

		if mod in self.models: return self.models[mod]
		return None

	def NewChain(self,name):

		""" Construct and add a new chain by name to the PDB. Returns the chain. """

		p = PDBchain(name)
		self.AddChain(name,p)
		return p

	def NewModel(self,name):

		""" Construct and add a new model by name to the PDB. Returns the model. """

		p = PDBmodel(name)
		self.AddModel(name,p)
		return p

	def RemoveChain(self,name):

		""" Remove a chain from the structure (by name). Returns the chain. """

		c = self.GetChain(name)
		if not name in self.chains:
			raise KeyError('Chain does not exist in structure by that name!')
		del self.chains[name]
		self.orderofchains.remove(name)
		return c

	def RemoveModel(self,name):

		""" Remove a model from the structure (by name). Returns the chain. """

		m = self.GetModel(name)
		if not name in self.models:
			raise KeyError('Model does not exist in structure by that name!')
		del self.models[name]
		self.orderofmodels.remove(name)
		return m

	def AddChain(self,chainname,chain):

		""" Add a chain as a list of residues to the PDB. """

		# Do chain addition operations.
		if chainname in self.chains:
			raise KeyError('Chain already exists in structure by that name!')
		if type(chain) != PDBchain and type(chain) == PDBmodel:
			# Cast into a PDBchain.
			cast = PDBchain(chainname)
			for i in chain: cast.AddResidue(chain[i])
			chain = cast
		self.chains[chainname] = chain
		self.orderofchains.append(chainname)
		chain.structure = self

	def AddModel(self,modelname,model):

		""" Add a model as a list of residues to the PDB. """

		# Do chain addition operations.
		if model in self.models:
			raise KeyError('Model already exists in structure by that name!')
		if type(model) != PDBmodel and type(model) == PDBchain:
			# Cast into a PDBmodel.
			cast = PDBmodel(modelname)
			for i in model: cast.AddResidue(model[i])
			model = cast
		self.models[modelname] = model
		self.orderofmodels.append(modelname)
		model.structure = self

	def AddResidueToChain(self, chain, res):

		""" Add a residue to a chain. Deprecated; use chain class function. """

		ch = self.GetChain(chain)
		if ch != None: ch.AddResidue(res)

	def AddResidueToModel(self, model, res):

		""" Add a residue to a model. Deprecated; use model class function. """

		mo = self.GetModel(model)
		if mo != None: mo.AddResidue(res)      

	def AddRemark(self,remark):

		""" Add a remark (note/comment) to the structure/PDB file. """

		if len(remark) == 0:
			self.remarks.append('')
			return
		for it in xrange(0,len(remark),68):
			self.remarks.append(remark[it:it+68])

	def GetRemarks(self):

		""" Return all remarks from the PDB as a list of strings. """

		return self.remarks

	def CheckComplete(self):

		""" For every chain, check to see if every residue is complete (see 
        aa_list dictionary). """

		for ch in self.orderofchains:
			for re in self.chains[ch]:
				res = self.chains[ch][re]
				nam = res.name
				if not nam in aa_lists: continue
				check = set(set(aa_lists[nam])).difference(res.atoms)
				if len(check) > 0:
					return res
		for ch in self.orderofmodels:
			for res in self.models[ch]:
				nam = res.name
				if not nam in aa_lists: continue
				check = set(set(aa_lists[nam])).difference(res.atoms)
				if len(check) > 0:
					return res                
		return True    

	# I/O Functionality.

	def _pymol(self):

		# Start pymol if not started.
		import pymol
		pymol.finish_launching()        
		return pymol

	def view(self,istrajectory=False):

		""" View the structure in a Pymol window. Requires an installation of Pymol. """

		pym = self._pymol()
		pym.cmd.read_pdbstr(str(self),'PDBnet')
		if len(self.models) > 0 and not istrajectory:
			pym.cmd.set('all_states','on')
		elif istrajectory: pass # SERGIO

	def ViewStructure(self): self.view()    

	def WriteFile(self, filename):

		""" Write this PDB structure as a single PDB file. """

		fh = open(filename,'w')
		fh.write(str(self))
		fh.close()

	def read(self,filename):

		""" Alias for ReadFile(). """

		self.ReadFile(filename)

	def ReadFile(self, filename):

		""" Read a PDB file. Populate this PDBstructure. """

		# Map to memory and read the file.
		p = PDBfile(filename)
		self.handle = p
		self.remarks, self.organism, self.taxid, self.mutation = p.read()
		largeFile = p.isLargeFile()

		# Start creating high-level structures for models/chains.
		if p.hasModels():
			self.ismodel = True
			map(self.NewModel,p.getModelNames())		
			for mo, ch, re in p.iterResidueData():
				model = self.GetModel(mo)
				if ch != '':
					if not ch in model.GetChainNames(): model.NewChain(ch)
					chain = model.GetChain(ch)
					if largeFile: chain.AddIndexOfResidue(re)
					else: chain.AddResidueByIndex(re)
				else:
					if largeFile: model.AddIndexOfResidue(re)
					else: model.AddResidueByIndex(re)
			for mod in self.GetModelNames():
				model = self.GetModel(mod)
				model.structure = self			
				model.SortByNumericalIndices()
				for chain in model.GetChains():
					chain.SortByNumericalIndices()		
		else:
			self.ismodel = False
			map(self.NewChain,p.getChainNames())			
			for _, ch, re in p.iterResidueData():
				chain = self.GetChain(ch)
				if largeFile: chain.AddIndexOfResidue(re)
				else: chain.AddResidueByIndex(re)
			for cha in self.GetChainNames():
				chain = self.GetChain(cha)
				chain.structure = self					
				chain.SortByNumericalIndices()

		# Remove this component from memory if no longer necessary.
		if not largeFile:
			p.close()
			self.handle = None
			del p

	def ChainAsFASTA(self, chain):

		""" Return the chain as a FASTA. Deprecated; use chain class function. """

		ch = self.GetChain(chain)
		if ch != None: return ch.AsFASTA()

	def ModelAsFASTA(self, model):

		""" Return the model as a FASTA. Deprecated; use chain class function. """

		mo = self.GetModel(model)
		if mo != None: return mo.AsFASTA()

	# Scoring Functionality

	def _iterResidueAssociations(self,fasta,chains=None,fast=None):

		""" PRIVATE: Use an input FASTA file to determine associations
        between residues as they exist between 2+ chains. Constructs a
		mapping between chain and sequence and yields all strictly
		homologous positions as residues, chain-by-chain."""

		ismodel=self.ismodel
		if ismodel:
			orderofthings = self.orderofmodels
			things = self.models
		else:
			orderofthings = self.orderofchains
			things = self.chains

		if not chains: chains = orderofthings
		chaininds = [orderofthings.index(x) for x in chains]  
		if len(chaininds) < 2:
			raise ValueError('Need more than two chains (gave %d).' % (len(chaininds)))

		if fast is None:
			fast = FASTAnet.FASTAstructure(fasta,uniqueOnly=False)
		seqs = fast.orderedSequences
		if len(things) != len(seqs):
			raise IOError('FASTA set length not equal to PDB length.')	
		posv = fast.getStrictlyUngappedPositions(chaininds)

		for n in chaininds: # Go chain by chain and yield information.
			seq = seqs[n].sequence
			ch  = orderofthings[n]
			matches = self.GetFASTAIndices(things[ch],seq)
			residueList = [match for match in matches]
			if residueList[0] is None:
				nm = seqs[n].name
				chn = things[ch].name
				raise IOError('No relation found between sequence %s and chain %s.' % (nm,chn))			
			residueList = sorted(residueList,key=lambda d:d.fstindex)
			for pos in posv: # Go position by position.
				yield pos, ch, residueList.pop(0) # Yield position, chain name, and residue.

	def tmscore(self,fasta,chains=None,native=None,CA=True):

		""" Get the TMscore between two chains. Requires a 
        FASTA alignment and a value for the length of the
        native structure (e.g., for a pairwise alignment,
        the length of the structure used as a reference
        before alignment was done). The latter is computed
        by assuming the first of both provided chains is the
        native structure; otherwise, uses a provided chain
        name (native input). """

		ismodel=self.ismodel
		if ismodel:
			orderofthings = self.orderofmodels
			things = self.models
		else:
			orderofthings = self.orderofchains
			things = self.chains

		if not chains: chains = orderofthings
		if len(chains) != 2:
			raise ValueError('Need exactly two chains to score.')
		fas     = FASTAnet.FASTAstructure(fasta,uniqueOnly=False)
		posvect = fas.getStrictlyUngappedPositions()
		items   = self._iterResidueAssociations(fasta,chains,fas)

		# Get the lengths of the respective chains.
		if native == None: leN = len(things[chains[0]]) 
		else: leN = len(things[native]) # Length of the reference structure.
		leT = len(posvect) # Number of aligned positions.

		# Calculate d_0 for this alignment.
		cuberoot = lambda d: d**(1./3)
		d_0 = 1.24 * cuberoot(leN-15) - 1.8

		# Get the summation portion of the TMscore.
		sumportion = 0
		posmask    = {}
		for pos,_,residue in items:
			if not pos in posmask: posmask[pos] = []
			if CA: atom = residue.GetCA()
			else:  atom = residue.Centroid()
			posmask[pos].append(atom)
		for pos in posvect: # Assume chA is Native Structure.
			cavect   = posmask[pos]
			ca1, ca2 = cavect[0], cavect[1]
			d_i      = sqrt((ca1.x-ca2.x)**2+
			                (ca1.y-ca2.y)**2+
			                (ca1.z-ca2.z)**2)
			sumdenom  = 1 + (d_i/d_0)**2
			suminside = 1./(sumdenom)
			sumportion += suminside

		# Return the TMscore.
		return (1./leN) * sumportion

	def gdt(self,fasta,chains=None,distcutoffs=[1,2,4,8],CA=True):

		""" Get the GDT score between two chains. Requires a
        FASTA alignment. """

		ismodel=self.ismodel

		if ismodel:
			orderofthings = self.orderofmodels
			things = self.models
		else:
			orderofthings = self.orderofchains
			things = self.chains

		if len(distcutoffs) != 4:
			raise ValueError('Need exactly 4 distance cutoff values.')
		if not chains: chains = orderofthings
		if len(chains) != 2:
			raise ValueError('Need exactly two chains to score.')
		items = self._iterResidueAssociations(fasta,chains)

		# Get all relevant atoms.
		posmask = {}
		for pos,_,residue in items:
			if not pos in posmask: posmask[pos] = []
			if CA: atom = residue.GetCA()
			else:  atom = residue.Centroid()
			posmask[pos].append(atom)		

		# Get the Euclidean distances between
		# all pairs of aligned positions.
		distances = {}
		for pos in posmask:
			cavect = posmask[pos]
			ca1, ca2 = cavect[0], cavect[1]
			dist = sqrt((ca1.x-ca2.x)**2+
			            (ca1.y-ca2.y)**2+
			            (ca1.z-ca2.z)**2)
			distances[pos] = dist

		# Calculate the counts for different cutoffs.
		counts  = [0,0,0,0]
		poslen  = float(len(distances))
		A,B,C,D = distcutoffs
		for pos in distances:
			dist = distances[pos]
			if dist < A: counts[0] += 1
			if dist < B: counts[1] += 1
			if dist < C: counts[2] += 1
			if dist < D: counts[3] += 1
		GDT_P1 = counts[0]/poslen
		GDT_P2 = counts[1]/poslen
		GDT_P3 = counts[2]/poslen
		GDT_P4 = counts[3]/poslen
		return (GDT_P1+GDT_P2+GDT_P3+GDT_P4)/4.0

	def rmsd(self,fasta,chains=None,CA=True):

		"""Get the RMSD between chains. Requires a FASTA alignment."""

		items = self._iterResidueAssociations(fasta,chains)
		# Get all relevant atoms.
		posmask = {}
		for pos,_,residue in items:
			if not pos in posmask: posmask[pos] = []
			if CA: atom = residue.GetCA()
			else:  atom = residue.Centroid()
			posmask[pos].append(atom)			
		# Get the Euclidean distance squared between
		# all pairs of aligned positions.
		distsqs = {}
		for pos in posmask:
			cavect = posmask[pos]
			for a in xrange(len(cavect)):
				for b in xrange(a+1,len(cavect)):
					if (a,b) not in distsqs: distsqs[(a,b)] = 0
					ca1, ca2 = cavect[a], cavect[b]
					distsq = ((ca1.x-ca2.x)**2
					          +(ca1.y-ca2.y)**2
					          +(ca1.z-ca2.z)**2)
					distsqs[(a,b)] += distsq
		# Calculate the average RMSD.
		rmsd = 0
		for d in distsqs:
			d = distsqs[d]
			r = sqrt((float(d)/len(posmask)))
			rmsd += r
		rmsd /= len(distsqs)
		return rmsd

	def rrmsd(self,fasta,chains=None,CA=True):

		""" Get the RRMSD between chains. Requires a FASTA alignment. 
        See Betancourt & Skolnick, "Universal Similarity Measure for
        Comparison Protein Structures". """

		if self.ismodel:
			orderofthings = self.orderofmodels
			things = self.models
		else:
			orderofthings = self.orderofchains
			things = self.chains

		if not chains: chains = orderofthings
		if len(chains) != 2:
			raise ValueError('Need exactly two chains to score.')        

		R_A = self.RadiusOfGyration([chains[0]])
		R_B = self.RadiusOfGyration([chains[1]])
		avglen = (len(things[chains[0]])+
		          len(things[chains[1]]))/2
		L_N, e = avglen-1, math.e
		c = 0.42-0.05*L_N*e**(-L_N/4.7)+0.63*e**(-L_N/37)
		denom = R_A**2+R_B**2-2*c*R_A*R_B
		alignRMSD = self.rmsd(fasta,chains,CA)
		return float(alignRMSD)/sqrt(denom)

	# Other Functionality

	def GetAverage(self,chains=None,newname=None):

		""" Acquire a new chain or model corresponding to an average of all present
        chains or models specified. """

		if self.ismodel:
			orderofthings = self.orderofmodels
			things        = self.models
			if newname == None: newname = 1
			output        = PDBmodel(newname)
		else:
			orderofthings = self.orderofchains
			things        = self.chains
			if newname == None: newname = '*'
			output        = PDBchain(newname)
		if chains == None: chains = orderofthings

		# Make sure all sequences have equal FASTAs.
		fas = ''
		for ch in chains:
			if fas == '': fas = things[ch].AsFASTA()
			if things[ch].AsFASTA() != fas:
				raise AssertionError('Not all provided chains/models have equal sequences.')

		def avgResidue(res_no):

			# Create the dummy residue structure.
			firstch  = things[chains[0]]
			firstres = firstch.GetResidues()[res_no]
			fake = PDBresidue(firstres.index,firstres.name)

			# Populate it with atoms.
			atomlist = []
			for a in firstres.GetAtoms():
				atomlist.append(PDBatom(a.serial,a.name,0,0,0,
				                        0,0,a.symbol,''))            

			# Average all properties.
			count = len(chains)
			for ch in chains:
				res = things[ch].GetResidueByIndex(things[ch].GetIndices()[res_no])
				atoms = res.GetAtoms()
				for i in xrange(len(atoms)):
					atomlist[i].x += atoms[i].x
					atomlist[i].y += atoms[i].y
					atomlist[i].z += atoms[i].z
					atomlist[i].occupancy += atoms[i].occupancy
					atomlist[i].tempFactor += atoms[i].tempFactor
			for a in atomlist:
				a.x /= float(count)
				a.y /= float(count)
				a.z /= float(count)
				a.occupancy /= float(count)
				a.tempFactor /= float(count)

			for a in atomlist: fake.AddAtom(a)
			return fake

		res_nums = range(len(things[chains[0]]))
		for i in res_nums: output.AddResidue(avgResidue(i))

		return output

	def RadiusOfGyration(self,chains=None):

		""" Acquire the radius of the gyration of the entire, or a portion of, the
        PDB protein molecule. """

		if self.ismodel:
			orderofthings = self.orderofmodels
			things = self.models
		else:
			orderofthings = self.orderofchains
			things = self.chains      
		if not chains: chains = orderofthings

		# Get centre of mass of entire molecule.
		molCentroid = [0.,0.,0.]
		numAtoms    = 0
		for c in chains:
			for re in things[c]:
				r = things[c][re]
				for at in r.atoms:
					a = r.atoms[at]
					molCentroid[0] += a.x
					molCentroid[1] += a.y
					molCentroid[2] += a.z
					numAtoms += 1
		m_x,m_y,m_z = molCentroid
		molCentroid = [m_x/numAtoms,m_y/numAtoms,m_z/numAtoms]

		# Sum all of the distances squared vs. centre of mass.
		distsqr = []
		m_x,m_y,m_z = molCentroid
		for c in chains:
			for re in things[c]:
				r = things[c][re]
				for at in r.atoms:
					a = r.atoms[at]
					diff = (a.x-m_x+a.y-m_y+a.z-m_z)
					distsqr.append(diff**2)

		# Return the sum of all of these, divided by number of atoms,
		# square rooted.
		sumdistratio = sum(distsqr)/float(numAtoms)
		return sqrt(sumdistratio)

	def GetAllCentroid(self, chain):

		""" Populates the centroids of all residues. """

		out = []
		if chain:
			if not chain in self.chains:
				raise ValueError('Chain %s does not exist!' % (chain))

			ch = self.chains[chain]

			for res in ch:
				res = ch[res]
				res.Centroid()
				out.append(res.centroid)
		else:
			for res in self.IterAllResidues():
				res.Centroid()
				out.append(res.centroid)				

		return out

	def IndexSeq(self, chain, fst):

		""" Store in residues the correct index to the fasta.
        Requires a 1-to-1 correspondance at least a portion
        of the way through. Deprecated; use GetFASTAIndices(). """

		ismodel=self.ismodel

		if ismodel: thing = self.GetModel(chain)
		else: thing = self.GetChain(chain)
		return self.GetFASTAIndices(thing, fst)            

	def GetFASTAIndices(self, thing, fst):

		""" Given a PDBchain, find 1-to-1 correspondances between
        it and a FASTA sequence. """

		chainseq = thing.AsFASTA()
		ungapped = fst.replace('-','')
		success  = True
		if len(ungapped) != len(chainseq):
			success = False
			yield None
		for i in xrange(0,len(chainseq)):
			# See if there is a failed correspondance.
			if chainseq[i] != ungapped[i]: 
				yield None
				success = False
				break
		if success:
			index = -1
			for i in thing:
				i = thing[i]
				index = fst.find(aa[i.name],index+1)
				i.fstindex = index
				yield i

	def IterResiduesFor(self,chains=None):

		""" Produce an iterator to allow one to iterate over all residues for a subset of the structure. """

		if self.ismodel:
			orderofthings = self.orderofmodels
			things        = self.models
		else:
			orderofthings = self.orderofchains
			things        = self.chains
		if chains == None: chains = orderofthings

		for ch in chains:
			chain = things[ch]
			for res in chain: yield chain[res]

	def IterAllResidues(self):

		""" Produce an iterator to allow one to iterate over all possible residues. """

		for ch in self.chains.keys():
			chain = self.GetChain(ch)
			for res in chain: yield chain[res]   
		for mo in self.models.keys():
			model = self.GetModel(mo)
			for res in model: yield model[res]

	def gm(self,fasta,chains=None,CA=False,typeof='str'):

		""" Acquire Geometric Morphometric data corresponding to the
        (x,y,z) coordinates between all homologous residue positions.
        Requires a FASTA alignment. Options include using alpha-carbon
        positions. By default, uses centroids of residues. Returns a list
        of labels and a list of coordinates as raw GM data. The typeof option 
        provides an option for coordinate output; they are returned 
        as a semicolon-delimited string (str) or as a numpy 2d array (matrix). """

		if self.ismodel:
			orderofthings = self.orderofmodels
			things        = self.models
		else:
			orderofthings = self.orderofchains
			things        = self.chains
		if chains == None: chains = orderofthings
		labels,coords = [],[] # Output Variables

		# Set up name grabbing from FASTA file for GM labelling.
		f = FASTAnet.FASTAstructure(fasta,uniqueOnly=False)
		name = lambda d: f.orderedSequences[d].name 

		# Acquire all homologous positions as defined in FASTA.
		items = self._iterResidueAssociations(fasta,chains,f)

		chainsSaw = []
		if typeof == 'matrix': row = []
		else:                  row = ''
		for pos,chain,res in items:
			if chain not in chainsSaw:
				ind = orderofthings.index(chain)
				nm = name(ind).strip().split(':')[0]
				labels.append('>%s:%s' % (nm,nm[:4]))
				if len(chainsSaw) != 0:
					coords.append(row)
				chainsSaw.append(chain)
				if typeof == 'matrix': row = []
				else: row = ''					
			if CA: atom = res.GetCA()
			else:  atom = res.Centroid()
			if isinstance(row, str):   row += '%f;%f;%f;'%(atom.x,atom.y,atom.z)
			elif isinstance(row,list): row.extend([atom.x,atom.y,atom.z])
		coords.append(row) # Last chain.

		if typeof == 'matrix': coords = np.array(coords)
		return labels,coords

	def WriteGM(self,fasta,gm,chains=None,CA=False):

		""" Write the information present in this PDB between multiple
        chains as a Geometric Morphometric text file. This file will be
        formatted such that individual lines correspond to chains and semi-colons
        separate the (x,y,z) coordinates between all homologous residue positions. 
        Requires a FASTA alignment. Options include using alpha-carbon positions.
        By default, uses centroids of residues. """

		fgm = open(gm,'w')

		# Get raw GM data.
		labels,coords = self.gm(fasta,chains,CA)

		# Start writing.
		for ind in xrange(len(labels)):
			fgm.write(labels[ind] + ';' + coords[ind] + '\n')

		# Done.
		fgm.close()

	def WriteLandmarks(self,fasta,lm,chains=None):

		""" Write the information present in this PDB between multiple
        chains as a landmark text file. This file will be formatted such that 
        the file is partitioned in sections starting with chain names and individual 
        lines in these correspond to homologous residue positions denoted by
        homologous position, residue number, and residue name tab-delimited. Requires
        a FASTA file. """

		ismodel=self.ismodel

		if ismodel:
			orderofthings = self.orderofmodels
			things        = self.models
		else:
			orderofthings = self.orderofchains
			things        = self.chains
		if chains == None: chains = orderofthings

		flm = open(lm,'w')

		# Set up name grabbing from FASTA file for labelling.
		f = FASTAnet.FASTAstructure(fasta,uniqueOnly=False)
		name = lambda d: f.orderedSequences[d].name

		# Start writing.
		chainsSaw = []
		ind       = -1
		for pos,chain,res in items:
			if chain not in chainsSaw:
				ind = orderofthings.index(chain)
				nm = name(ind).strip().split(':')[0]
				flm.write('>%s\n' % (nm))
				chainsSaw.append(chain)
				ind = 0
			flm.write('%s\t%s\t%s\n' % (ind,res.index,res.name))
			ind += 1

	def FDmatrix(self,fasta,chains=None,scaled=True):
		"""
		Compute the form difference matrix (FDM) as explained in Claude 2008. It relates to
		the identification of the most influential residue, with respect to the overall
		shape/structure. If the scaled option is True, will return an scaled list 
		(from -1 to 1) of the of lenght equal to the number of residues. Otherwise will return
		the raw FDM, rounded so it can be included in a PDB. The scaled version is better for
		vizualization. By default the FDM is computed far all chains, but a subset can be passed
		to the chains option.
		"""
		names, data = self.gm(fasta,chains=chains,typeof='matrix')
		#get dimensions
		n_obs = data.shape[0]     # number observations/structures
		n_res = data.shape[1]/3   # Number of residues		
		# Re-sahpe data as a 3D array
		data= data.reshape(n_obs,n_res,3)
		# compute meanshape	
		msh=data.mean(axis=0)
		#compute the form matrix for the mean shape
		FMmsh = cdist(msh,msh)
		FMmsh = np.ma.masked_array(FMmsh,np.isnan(FMmsh))
		#iterate over observations to compute the form difference matrix
		vector=[]
		for index in xrange(n_obs):
			ndarray=data[index,:,:]
			FM   = cdist(ndarray,ndarray)
			FDM  = np.ma.masked_array(FM,np.isnan(FM))/FMmsh
			utr  = FDM[np.triu_indices(FDM.shape[0],k=-1)]
			mfdm = np.absolute(FDM - np.median(utr[~np.isnan(utr)]))
			mdat = np.ma.masked_array(mfdm,np.isnan(mfdm))
			s    = np.sum(mdat,axis=1).filled(np.nan)
			vector.append(s)

		#compute FD and scaled if required
		FD  = (np.array(vector).sum(axis=0))/n_obs
		if not scaled:
			return [round(x,2) for x in (FD)]
		else:
			return [round(x,3) for x in (2*(FD - min(FD))/(max(FD) - min(FD)) - 1)]

	def Map2Protein(self,outname,lis,chain,fasta):
		"""
		Map a list of values (lis), that must have a lenght equal to that of the number
		of residues in the PDB to be mapped (chain). If a list of list is provided, the 
		first list will be mapped as the beta factor and the second as occupancy
		"""
		# Check if more than one thing need to be included
		if len([lis]) == 1:
			lis=[lis]

		# Construct empty PDBstructure.
		dummy = PDBstructure()

		# Reset beta and ocuppancy
		if self.models:
			m = self.GetModel(chain)
			newm = dummy.NewModel(m.name)
			for r in m.IterResidues():
				for a in r.GetAtoms():
					a.tempFactor=0.00
					if len(lis) > 1:
						a.occupancy=0.00
				newm.AddResidue(r)
		else:
			m = self.GetChain(chain)
			newm = dummy.NewChain(m.GetName())
			for r in m.IterResidues():
				for a in r.GetAtoms():
					a.tempFactor=0.00
					if len(lis) > 1:
						a.occupancy=0.00
				newm.AddResidue(r)			

		# Acquire all homologous positions as defined in FASTA.
		items = self._iterResidueAssociations(fasta,None,None)		

		# populate new field
		for pos,chainb,res in items:
			if not chainb == chain:
				continue
			res = newm.GetResidueByIndex(res.index)
			for a in res.GetAtoms():
				a.tempFactor = lis[0][pos]
				if len(lis) > 1:
					a.occupancy = lis[1][pos]


		newm.WriteAsPDB(outname)


	def Contacts(self,chain=None,thres=4.5):

		""" Compute the contact map of all chains or a chain.

        :param chain: A list of chain or model names or a single string or integer. By default, entire structure.
        :param thres: A threshold for distinguishing contact in Angstroms.
        :returns: A list of tuples of indices (integers) which correspond to chains or models and their residues. 

        """     

		ismodel = self.ismodel

		if ismodel:
			orderofthings = self.orderofmodels
			things        = self.models
		else:
			orderofthings = self.orderofchains
			things        = self.chains

		if chain == None or chain == 'ALL':
			chain = orderofthings

		done            = []
		self.contactmap = []
		if (type(chain) == int or type(chain) == str) or len(chain) == 1:
			if type(chain) == list or type(chain) == tuple:
				chain = things[chain[0]]
			else: chain = things[chain]
			self.contactmap = chain.ContactMap(thres)
		else: # A number of chains.
			iteratorA = self.IterResiduesFor(chain)
			iteratorB = self.IterResiduesFor(chain)
			for resA in iteratorA:
				for resB in iteratorB:
					if resA == resB or (resA,resB) in done: continue 
					elif resA.InContactWith(resB,thres):
						if not ismodel:
							A = resA.index+resA.chain
							B = resB.index+resB.chain
						else:
							A = (resA.index,resA.model)
							B = (resB.index,resB.model)
						if not (A,B) in self.contactmap:
							self.contactmap.append((A,B))
							self.contactmap.append((B,A))
					done.append((resA,resB))
					done.append((resB,resA))
			self.contactmap = sorted(self.contactmap, key=lambda element: (
			    element[0], element[1]))

		return self.contactmap

	def WriteContacts(self, filename):

		""" Write contact map. """

		fout = open(filename, 'w')
		for a in self.contactmap:
			fout.write('%s\n'%(str(a)))
		fout.close()            

	# Internals

	def __iter__(self):

		""" Returns all PDBchains as an iterator. """

		for ch in self.chains: yield self.chains[ch]
		for mo in self.models: yield self.models[mo]

	def __len__(self):

		""" Returns the length of the protein. """

		# For all chains.
		chs = self.chains
		return sum([len(chs[x]) for x in chs])

	def _remarksToString(self):

		remarkstr = ''
		for it in xrange(len(self.remarks)):
			remarkstr += 'REMARK%4s %s\n' % (str(it+1),self.remarks[it])
		return remarkstr

	def __str__(self):

		""" As a string, outputs structure in the PDB format. """

		out = ''
		out += self._remarksToString()
		for chain in self.orderofchains: out += str(self.GetChain(chain))
		for model in self.orderofmodels: out += str(self.GetModel(model))
		return out

class PDBfile(object):

	__slots__ = ('filePath','fileHandle','memHandle','modelIndices','chains','residueIndices','size')

	def __init__(self,fi):

		self.filePath       = fi
		self.fileHandle     = open(fi,'r+b')
		self.size           = sizeOfFile(fi)
		self.memHandle      = memoryMappedFile(self.fileHandle.fileno(),0)
		self.residueIndices = dict()
		self.modelIndices   = dict()
		self.chains         = list()

	def isLargeFile(self):

		""" Return whether or on the PDB file this object represents is incredibly large. """

		return (self.size >= PDB_LARGE_FILE_SIZE)	

	def close(self):

		""" Close the file. """

		if self.memHandle:
			self.memHandle.close()
		else: self.fileHandle.close()

	def read(self):

		""" Acquire all remarks, and the indices of all models and residues. Returns
		remarks, and biological source information as a tuple (remarks, organism, 
		taxid, mutant). """

		# Metadata.
		remarks,organism,taxid,mutant = [],'','',False

		# Cursors.
		k      = 0
		curRes = None
		pos    = 0
		line   = self.memHandle.readline()
		model  = -1

		# Scan the file.
		while line:
			recordType = line[:6].strip()
			if recordType == 'MODEL':
				model = int(line.strip().split()[-1])
				self.modelIndices[model] = pos
			elif recordType == 'ATOM':
				sl = line.split()
				if int(sl[1]) >= 100000: k = 1
				else: k = 0							
				chain = line[21+k].strip()
				if chain not in self.chains:
					self.chains.append(chain)
				residue_index = line[22+k:27+k].strip()
				iCode = line[27].strip()
				residue = residue_index + iCode
				if (curRes == None or residue != curRes):
					self.residueIndices[(model,chain,residue)] = pos
					curRes = residue
			elif recordType == 'REMARK':
				remark = line[6:].strip('\n')
				if not remark in remarks: remarks.append(remark)
			elif recordType == 'SOURCE':
				if 'ORGANISM_SCIENTIFIC' in line:
					organism = '_'.join(line.strip().strip(';').split()[-2:])
				elif 'ORGANISM_TAXID' in line:
					taxid    = line.strip().strip(';').split()[-1]
				elif 'MUTANT' in line or 'MUTATION' in line:
					mutant   = True
			pos  = self.memHandle.tell()
			line = self.memHandle.readline()
		self.memHandle.seek(0)

		# Return metadata.
		return (remarks,organism,taxid,mutant)

	def readResidue(self,chain,res,model=-1):

		""" Parse a residue from the PDB file and return a PDBresidue. """

		resObj = None
		index = self.residueIndices[(int(model),chain,res)]
		self.memHandle.seek(index)
		line = self.memHandle.readline()
		while line:
			recordType = line[:6].strip()
			if recordType == 'ATOM':
				sl = line.split()
				if int(sl[1]) >= 100000: k = 1
				else: k = 0		
				chain      = line[21+k].strip()
				residue    = line[22+k:27+k].strip() + line[27].strip()
				residue_id = line[17+k:20+k].strip()
				if residue == res:
					if not resObj:
						resObj = PDBresidue(residue,residue_id)
					serial = int(line[6+k:12+k])
					atom_name = line[12+k:16+k].strip()
					x = float(line[30+k:38+k])
					y = float(line[38+k:46+k])
					z = float(line[46+k:54+k])
					o = float(line[54+k:60+k])
					b = float(line[60+k:66+k])
					sym = line[76+k:78+k].strip()
					try: charge = line[79+k:80+k].strip()
					except: charge=' '					
					atom = PDBatom(serial,atom_name,x,y,z,o,b,sym,charge)
					resObj.AddAtom(atom)
				else: break
			line = self.memHandle.readline()
		self.memHandle.seek(0)
		return resObj

	def hasModels(self): return (len(self.modelIndices) > 0)
	def getModelNames(self): return self.modelIndices.keys()
	def getChainNames(self): return self.chains

	def iterResidueData(self):

		""" Yield the model number, chain name, and residue number 
		for all residues present in this PDB file, not necessarily 
		in order. """

		for mod,ch,res in self.residueIndices.keys():
			yield mod,ch,res

	def getResidueNamesInChain(self,ch):
		res = []
		for c,r in self.residueIndices.keys():
			if c == ch: res.append(r)
		return res

# Constants

aa = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E',
      'PHE':'F','GLY':'G','HIS':'H','ILE':'I',
      'LYS':'K','LEU':'L','MET':'M','ASN':'N',
      'PRO':'P','GLN':'Q','ARG':'R','SER':'S',
      'THR':'T','VAL':'V','TRP':'W','TYR':'Y',
      'UNK':'X'}

aa_names = {v:k for k, v in aa.items()} # Flip above dictionary.

aa_lists = {'ALA':['N','CA','C','O','CB'],\
            'CYS':['N','CA','C','O','CB','SG'],\
            'ASP':['N','CA','C','O','CB','CG','OD1','OD2'],\
            'GLU':['N','CA','C','O','CB','CG','CD','OE1','OE2'],\
            'PHE':['N','CA','C','O','CB','CG','CD1','CD2','CE1','CE2','CZ'],\
            'GLY':['N','CA','C','O'],\
            'HIS':['N','CA','C','O','CB','CG','ND1','CD2','CE1','NE2'],\
            'ILE':['N','CA','C','O','CB','CG1','CG2','CD1'],\
            'LYS':['N','CA','C','O','CB','CG','CD','CE','NZ'],\
            'LEU':['N','CA','C','O','CB','CG','CD1','CD2'],\
            'MET':['N','CA','C','O','CB','CG','SD','CE'],\
            'ASN':['N','CA','C','O','CB','CG','OD1','ND2'],\
            'PRO':['N','CA','C','O','CB','CG','CD'],\
            'GLN':['N','CA','C','O','CB','CG','CD','OE1','NE2'],\
            'ARG':['N','CA','C','O','CB','CG','CD','NE','CZ','NH1','NH2'],\
            'SER':['N','CA','C','O','CB','OG'],\
            'THR':['N','CA','C','O','CB','OG1','CG2'],\
            'VAL':['N','CA','C','O','CB','CG1','CG2'],\
            'TRP':['N','CA','C','O','CB','CG','CD1','CD2','NE1','CE2','CE3','CZ2','CZ3','CH2'],\
            'TYR':['N','CA','C','O','CB','CG','CD1','CD2','CE1','CE2','CZ','OH'],\
            'UNK':[]}

# Debugging

if __name__ == "__main__":
	mystructure = PDBstructure(sys.argv[1])
	if mystructure.ismodel:
		names = mystructure.GetModelNames()
	else:
		names = mystructure.GetChainNames()
	lis = mystructure.FDmatrix(sys.argv[1].strip('pdb')+'gm') # for debugging purposes only
	mystructure.Map2Protein('test.pdb',lis,names[0],sys.argv[1])# for debugging purposes only
	if len(sys.argv) > 2:
		print 'RMSD',mystructure.rmsd(sys.argv[2])
		print 'RRMSD',mystructure.rrmsd(sys.argv[2])
		print 'TMscore',mystructure.tmscore(sys.argv[2])
		print 'GDT',mystructure.gdt(sys.argv[2])
	x =  mystructure.GetAllCentroid('A')
	mystructure.Contacts()
	print mystructure.contactmap