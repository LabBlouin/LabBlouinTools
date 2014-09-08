#!/bin/python

''' 
PDBnet is a collection of Python objects intended to model and contain
PDB protein data.

PDBnet Copyright (C) 2012 Christian Blouin
Later Revisions: Christian Blouin, Alex Safatli and Jose Sergio Hleap

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

E-mail: cblouin@cs.dal.ca
Dependencies: BioPython, FASTAnet (contained in LabBlouinTools)

'''

import sys, FASTAnet, math
from math import sqrt
from Bio.SeqUtils.ProtParam import ProteinAnalysis as PA

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

# Classes

class PDBatom(object):
    
    ''' ATOM in a PDB protein structure. '''
    
    def __init__(self, serial, name, x,y,z, oc, b, symbol,charge):
        
        ''' Construct an atom in a PDB protein structure. '''
        
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
        
        ''' Ensures the name of this atom fits within 4 characters. '''
        
        if len(self.name) == 4: return self.name
        else: return ' '+self.name.ljust(3)
        
    def __str__(self):
        
        ''' Express this atom as a single line in a PDB file (see PDB format). '''
        
        return "ATOM  %5s %4s %3s %1s%4s    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s" \
               % (self.serial,self.fixname(),self.parent.name,self.parent.chain,
                  self.parent.index,self.x,self.y,self.z,self.occupancy,
                  self.tempFactor,self.symbol,self.charge)

    def GetPosition(self): 
        
        ''' Get the 3-dimensional coordinate of this atom in space. '''
        
        return (self.x,self.y,self.z)

    def DistanceTo(self, atom):

        ''' Acquire the distance from this atom to another atom. '''
        
        return ((self.x-atom.x)**2+(self.y-atom.y)**2+(self.z-atom.z)**2)**0.5

class PDBterminator(PDBatom):
    
    ''' A placeholder class that represents a terminating ATOM-like line in the PDB file. '''
    
    def __init__(self,chaininst):
        
        self.lastatom = chaininst.GetAtoms()[-1]
        self.lastresname = self.lastatom.parent.name
        self.lastreschain =  self.lastatom.parent.chain
        self.lastresind = self.lastatom.parent.index
        newserial = str(int(self.lastatom.serial) + 1)
        super(PDBterminator,self).__init__(newserial,'',0,0,0,0,0,'',0)

    def __str__(self):
        
        return 'TER   %5s %4s %3s %1s%4s' % (
            self.serial,'',self.lastresname,self.lastreschain,self.lastresind)

class PDBresidue:
    
    ''' A residue (collection of ATOM fields) in a PDB protein structure. '''
    
    def __init__(self,index=None,name=''):
        
        ''' Construct a residue (collection of atoms) in a PDB protein structure. '''
        
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
        
        ''' Add a PDBatom structure to this residue. '''
        
        self.atoms[atom.name.strip()] = atom
        self.atomsOrdered.append(atom)
        atom.parent = self

    def GetCA(self):
        
        ''' Get the alpha-carbon found in this residue as a PDBatom. '''
        
        if 'CA' in  self.atoms: return self.atoms['CA']
        else: return None

    def Centroid(self):
        
        ''' Calculate the centroid of this residue. Return this as a PDBatom. '''
        
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
        
        ''' Determine if in contact with another residue. '''
        
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
    
    ''' A PDB chain (collection of protein residues). '''
    
    def __init__(self,name):
        
        ''' Construct a PDB chain (collection of protein residues). '''
        
        self.name       = name
        self.structure  = None
        self.residues   = []
        self.indices    = []
    
    def GetResidues(self): return [self.GetResidueByIndex(i) for i in self.GetIndices()]
    def GetResidueByIndex(self,index): return self.__getitem__(index)
    
    def GetAtoms(self): 
        
        ''' Get all comprising atoms in this chain in order of residues. '''
        
        atoms = []
        for x in self.GetResidues(): atoms.extend(x.GetAtoms())
        return atoms
    
    def GetIndices(self): return [int(x) for x in self.indices]    
    
    def AddResidue(self,resid):

        ''' Add a residue to this chain. '''
        
        resid.chain = self.name
        self.residues.append(resid)
        self.indices.append(str(resid.index))
        
    def RemoveResidue(self,resid):
        
        ''' Remove a residue from this chain. '''
        
        return self.pop(resid)
    
    def GetPrimaryPropertiesFromBioPython(self):
        
        ''' Use BioPython to populate this class with attributes. '''
        
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
        
        ''' Return the string representing this chain's FASTA sequence. '''
        
        out = ''
        for item in self.residues: out += aa[item.name]
        return out

    def WriteAsPDB(self,filename):
        
        ''' Write this single chain to a file as a PDB. '''
        
        fh = open(filename,'w')
        fh.write('\n'.join([str(x) for x in self.residues]))
        fh.close()        
        
    def SortByNumericalIndices(self):
        
        ''' Sort all internal items by a numerical index. '''
        
        self.residues = sorted(self.residues,key=lambda d: int(d.index))
        self.indices  = sorted(self.indices, key=lambda d: int(d))

    def pop(self,s):
        
        ''' Mimic dictionary class function to pop a residue from this chain. '''
        
        res = self.GetResidueByIndex(s)
        self.residues.remove(res)
        self.indices.remove(str(s))
        return res

    def update(self,o):
        
        ''' Mimic dictionary class function to add residues from another chain. '''
        
        for res in o: self.AddResidue(res)

    def __str__(self):
        
        get = lambda d: self.GetResidueByIndex(d)
        out = ''
        out += '\n'.join([str(get(x)) for x in self.indices])
        out += '\n%s\n' % (str(PDBterminator(self)))
        return out

    def __getitem__(self,i):
        if str(i) in self.indices:
            return self.residues[self.indices.index(str(i))]
        else: return None

    def __len__(self): return len(self.residues)

    def __iter__(self):
        for ind in self.indices: yield int(ind)

class PDBmodel(PDBchain):
    
    ''' A PDB model (a special kind of chain). '''
    
    def __init__(self,name):
        
        ''' Construct a PDB model (a specific kind of chain). '''
        
        super(PDBmodel,self).__init__(name)
        self.chains = []

    def GetChain(self,name):
        return self.GetChainByName(name)
        
    def GetChainByName(self,name):
        for x in self.chains:
            if x.name == name: return x

    def GetChainNames(self): return [x.name for x in self.chains]

    def AddChain(self,chain):
        
        ''' Add a PDBchain (chain) to this model. '''

        self.chains.append(chain)

    def AddResidue(self,resid):
        
        ''' Add a residue to this model. '''
        
        resid.model = self.name
        self.residues.append(resid)
        self.indices.append(str(resid.index))
        
    def __str__(self):
        
        return 'MODEL%9s\n' % (str(self.name)) + super(
            PDBmodel,self).__str__() + 'ENDMDL\n'

class PDBstructure:
    
    ''' A PDB protein structure. Allows one to read in and parse PDB files. '''
    
    def __init__(self, filein=''):
    
        ''' Construct a protein PDB (a collection of chains). '''
        
        # Attributes.
        self.chains        = {}
        self.orderofchains = []
        self.models        = {}
        self.orderofmodels = []
        self.remarks       = []
        self.filepath      = filein
        self.organism      = None
        self.taxid         = None
        self.mutation      = False
        self.contactmap    = list()   
        
        # Read a file?
        if filein != '':
            self.ReadFile(filein)
        
    # Accesssors/Mutators.
    
    def GetChainNames(self): return self.orderofchains
    def GetModelNames(self): return self.orderofmodels
    
    def GetChain(self,ch): 
        
        ''' Get a chain by name. '''
        
        if ch in self.chains: return self.chains[ch]
        return None
    
    def GetModel(self,mod):
        
        ''' Get a model by name. '''
        
        if mod in self.models: return self.models[mod]
        return None

    def NewChain(self,name):
        
        ''' Construct and add a new chain by name to the PDB. Returns the chain. '''
        
        p = PDBchain(name)
        self.AddChainToStructure(name,p)
        return p
    
    def NewModel(self,name):
        
        ''' Construct and add a new model by name to the PDB. Returns the model. '''
        
        p = PDBmodel(name)
        self.AddModelToStructure(name,p)
        return p

    def AddChainToStructure(self,chainname,chain):

        ''' Add a chain as a list of residues to the PDB. '''

        # Do chain addition operations.
        if chainname in self.chains:
            raise AssertionError('Chain already exists in structure by that name!')
        if type(chain) != PDBchain and type(chain) == PDBmodel:
            # Cast into a PDBchain.
            cast = PDBchain(chainname)
            for i in chain: cast.AddResidue(chain[i])
            chain = cast
        self.chains[chainname] = chain
        self.orderofchains.append(chainname)
        chain.structure = self

    def AddModelToStructure(self,modelname,model):

        ''' Add a model as a list of residues to the PDB. '''

        # Do chain addition operations.
        if model in self.models:
            raise AssertionError('Model already exists in structure by that name!')
        if type(model) != PDBmodel and type(model) == PDBchain:
            # Cast into a PDBmodel.
            cast = PDBmodel(modelname)
            for i in model: cast.AddResidue(model[i])
            model = cast
        self.models[modelname] = model
        self.orderofmodels.append(modelname)
        model.structure = self

    def AddResidueToChain(self, chain, res):
        
        ''' Add a residue to a chain. Deprecated; use chain class function. '''
        
        ch = self.GetChain(chain)
        if ch != None: ch.AddResidue(res)

    def AddResidueToModel(self, model, res):
        
        ''' Add a residue to a model. Deprecated; use model class function. '''
        
        mo = self.GetModel(model)
        if mo != None: mo.AddResidue(res)      
    
    def AddRemark(self,remark):
        
        ''' Add a remark (note/comment) to the structure/PDB file. '''
        
        if len(remark) == 0:
            self.remarks.append('')
            return
        for it in xrange(0,len(remark),68):
            self.remarks.append(remark[it:it+68])
            
    def GetRemarks(self): return self.remarks
    
    def CheckComplete(self):
        
        ''' For every chain, check to see if every residue is complete (see 
        aa_list dictionary). '''
        
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
        
        ''' View the structure in a Pymol window. Requires an installation of Pymol. '''
        
        pym = self._pymol()
        pym.cmd.read_pdbstr(str(self),'PDBnet')
        if len(self.models) > 0 and not istrajectory:
            pym.cmd.set('all_states','on')
        elif istrajectory: pass # SERGIO
    
    def ViewStructure(self): self.view()    
        
    def write(self,filename):
        
        ''' Alias to WriteFile(). '''
        
        self.WriteFile(filename)    

    def WriteChain(self, ch, filename):
        
        ''' Write a single chain to a file. Deprecated; use chain class function. '''
        
        self.chains[ch].WriteAsPDB(filename)

    def WriteModel(self, mod, filename):
        
        ''' Write a single model to a file. Deprecated; use model class function. '''
        
        self.models[mod].WriteAsPDB(filename)

    def WriteFile(self, filename):
        
        ''' Write this PDB structure as a single PDB file. '''
        
        fh = open(filename,'w')
        fh.write(str(self))
        fh.close()

    def read(self,filename):
        
        ''' Alias for ReadFile(). '''
        
        self.ReadFile(filename)

    def ReadFile(self, filename):
        
        ''' Read a PDB file. Populate this PDBstructure. '''

        # Warning flag.
        hadwarning = False

        # Dummy PDB residue
        currentRes = PDBresidue()

        with open(filename) as fin:
        
            # Keep track of models (special type of chain).
            isModel = False
            model   = 0
        
            for line in fin:
                
                # Acquire organism and mutation metadata.
                if 'ORGANISM_SCIENTIFIC:' in line and not self.organism:
                    self.organism = '_'.join(line.strip().strip(';').split()[-2:])
                elif 'ORGANISM_TAXID:' in line and not self.taxid:
                    self.taxid = line.strip().strip(';').split()[-1]
                elif 'MUTANT' in line or 'MUTATION' in line:
                    self.mutation = True
                    
                # Store any appropriate REMARK entries.
                if line.startswith('REMARK'):
                    self.remarks.append(line[6:].strip('\n'))
                    continue
                    
                # Skip non atom/model entries in the file
                elif line.startswith('MODEL'):
                    model = int(line.strip().split()[-1])
                    isModel = True
                    continue
                elif line.startswith('ENDMDL'):
                    isModel = False
                    continue
                elif not line.startswith('ATOM'): continue

                # Break the line into fields.
                sl = line.split()
                if int(sl[1]) >= 100000: k = 1
                else: k = 0

                # Extract the useful information.
                chain = line[21+k].strip()
                serial = int(line[6+k:12+k])#modified by JSH
                residue_index = line[22+k:27+k].strip()#modified by JSH
                residue_id = line[17+k:20+k].strip()
                if residue_id == 'UNK' and not hadwarning:
                    hadwarning = True
                atom_name = line[12+k:16+k].strip()
                alt_loc=line[16].strip()
                iCode = line[27].strip()
                x = float(line[30+k:38+k])
                y = float(line[38+k:46+k])
                z = float(line[46+k:54+k])
                occupancy = float(line[54+k:60+k])#modified by JSH
                temp = float(line[60+k:66+k])#modified by JSH
                symbol = line[76+k:78+k].strip()#modified by JSH
                try: charge = line[79+k:80+k].strip()#modified by JSH
                except: charge=' '

                # Ensure inclusion of iCode in residue_index.
                residue = residue_index + iCode

                # populate model if MD
                if isModel:
                    if not model in self.models: self.NewModel(model)
                    if chain != '' and not chain in self.GetModel(model).GetChainNames():
                        # Is a MODEL and also has a chain.
                        self.GetModel(model).AddChain(PDBchain(chain))

                # Create chain instance if needed
                else:
                    if not chain in self.chains: self.NewChain(chain)

                # Create residue instance if needed
                if currentRes.index != residue:
                    currentRes = PDBresidue(residue,residue_id)
                    if not isModel:self.AddResidueToChain(chain,currentRes)
                    else: 
                        if chain != '':
                            # Is a MODEL and also has a chain.
                            self.GetModel(model).GetChainByName(chain).AddResidue(currentRes)
                        self.AddResidueToModel(model,currentRes)

                # Add atom to currentRes
                currentRes.AddAtom(PDBatom(
                    serial,atom_name,x,y,z,occupancy,temp,symbol,charge))

        # Return if had warnings.
        return hadwarning

    def ChainAsFASTA(self, chain):
        
        ''' Return the chain as a FASTA. Deprecated; use chain class function. '''

        ch = self.GetChain(chain)
        if ch != None: return ch.AsFASTA()

    def ModelAsFASTA(self, model):

        ''' Return the model as a FASTA. Deprecated; use chain class function. '''

        mo = self.GetModel(model)
        if mo != None: return mo.AsFASTA()

    # Scoring Functionality

    def GetResidueAssociations(self,fasta,chains=None,ismodel=False):

        ''' Use an input FASTA file to determine associations
        between residues as they exist between 2+ chains. Returns
        a mapping between chain and sequence, a masking of gapped
        parts of the sequence(s), and a list of aligned positions.'''
        
        if ismodel:
            orderofthings = self.orderofmodels
            things = self.models
        else:
            orderofthings = self.orderofchains
            things = self.chains

        if not chains: chains = orderofthings
        chaininds = [orderofthings.index(x) for x in chains]  
        if len(chaininds) < 2:
            raise ValueError('Need more than two chains.')
            
        fast = FASTAnet.FASTAstructure(fasta,uniqueOnly=False)
        seqs = fast.orderedSequences
        if len(things) != len(seqs):
            raise IOError('FASTA set length not equal to PDB length.')
        
        # Build a set of masks to represent the FASTA
        # sequence alignment as corresponding to the
        # PDB.
        mapng = {} # Mapping between chain name and sequence index.
        masks = {} # Mapping between residue and position.
        for n in chaininds:
            seq = seqs[n].sequence
            ch  = orderofthings[n]
            match = self.GetFASTAIndices(things[ch],seq)
            if not match: raise IOError(
                'Seq %d could not map to the respective PDB chain' % (n))
            else: mapng[ch] = n
            relist = [things[ch][x] for x in things[ch]]
            relist = sorted(relist,key=lambda d:d.fstindex)
            masks[n] = []
            for i in xrange(len(seq)):
                char = seq[i]
                if char.isalpha():
                    masks[n].append(relist.pop(0))
                else: masks[n].append(None)
                
        # Go through masks and determine what
        # positions are fully aligned.
        posvect = []
        alnlen  = len(masks.values()[0])
        for pos in xrange(alnlen):
            all_aligned = True
            for chain in mapng:
                if not masks[mapng[chain]][pos]:
                    all_aligned = False
                    break
            if all_aligned: posvect.append(pos) 

        return mapng, masks, posvect

    def tmscore(self,fasta,chains=None,ismodel=False,native=None):

        ''' Get the TMscore between two chains. Requires a 
        FASTA alignment and a value for the length of the
        native structure (e.g., for a pairwise alignment,
        the length of the structure used as a reference
        before alignment was done). The latter is computed
        by assuming the first of both provided chains is the
        native structure; otherwise, uses a provided chain
        name (native input). '''

        if ismodel:
            orderofthings = self.orderofmodels
            things = self.models
        else:
            orderofthings = self.orderofchains
            things = self.chains
            
        if not chains: chains = orderofthings
        if len(chains) != 2:
            raise ValueError('Need exactly two chains to score.')
        mapng,masks,posvect = self.GetResidueAssociations(fasta,chains,ismodel)

        # Get the lengths of the respective chains.
        chA,chB = mapng # Chains associated with alignment.
        if native == None: leN = len(things[chA]) 
        else: leN = len(things[native]) # Length of the reference structure.
        leT     = len(posvect) # Number of aligned positions.

        # Calculate d_0 for this alignment.
        cuberoot = lambda d: d**(1./3)
        d_0 = 1.24 * cuberoot(leN-15) - 1.8

        # Get the summation portion of the TMscore.
        sumportion = 0
        for pos in posvect: # Assume chA is Native Structure.
            cavect = []
            for ch in chA,chB:
                cavect.append(masks[mapng[ch]][pos].GetCA())
            ca1, ca2 = cavect[0], cavect[1]
            d_i       = sqrt((ca1.x-ca2.x)**2+
                             (ca1.y-ca2.y)**2+
                             (ca1.z-ca2.z)**2)
            sumdenom  = 1 + (d_i/d_0)**2
            suminside = 1./(sumdenom)
            sumportion += suminside

        # Return the TMscore.
        return (1./leN) * sumportion

    def gdt(self,fasta,chains=None,distcutoffs=[1,2,4,8],ismodel=False):

        ''' Get the GDT score between two chains. Requires a
        FASTA alignment. '''
        
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
        mapng,masks,posvect = self.GetResidueAssociations(fasta,chains,ismodel)

        # Get the Euclidean distances between
        # all pairs of aligned positions.
        distances = {}
        for pos in posvect:
            cavect = []
            for ch in mapng: 
                cavect.append(masks[mapng[ch]][pos].GetCA())
            if len(cavect) != 2: raise AssertionError(
                'Too many CA vectors.')
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

    def rmsd(self,fasta,chains=None,ismodel=False):

        '''Get the RMSD between chains. Requires a FASTA alignment.'''

        mapng,masks,posvect = self.GetResidueAssociations(fasta,chains,ismodel)
        # Get the Euclidean distance squared between
        # all pairs of aligned positions.
        distsqs = {}
        for pos in posvect:
            cavect = []
            for ch in mapng: 
                cavect.append(masks[mapng[ch]][pos].GetCA())
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
            r = sqrt((float(d)/len(posvect)))
            rmsd += r
        rmsd /= len(distsqs)
        return rmsd

    def rrmsd(self,fasta,chains=None,ismodel=False):
        
        ''' Get the RRMSD between chains. Requires a FASTA alignment. 
        See Betancourt & Skolnick, "Universal Similarity Measure for
        Comparison Protein Structures". '''

        if ismodel:
            orderofthings = self.orderofmodels
            things = self.models
        else:
            orderofthings = self.orderofchains
            things = self.chains
        
        if not chains: chains = orderofthings
        if len(chains) != 2:
            raise ValueError('Need exactly two chains to score.')        
        
        def getAverageCorrelationCoefficient(len_struct):
            L_N, e = len_struct-1, math.e
            return 0.42-0.05*L_N*e**(-L_N/4.7)+0.63*e**(-L_N/37)        
        
        R_A = self.RadiusOfGyration([chains[0]],ismodel)
        R_B = self.RadiusOfGyration([chains[1]],ismodel)
        avglen = (len(things[chains[0]])+
                  len(things[chains[1]]))/2
        c = getAverageCorrelationCoefficient(avglen)
        denom = R_A**2+R_B**2-2*c*R_A*R_B
        alignRMSD = self.rmsd(fasta,chains,ismodel)
        return float(alignRMSD)/sqrt(denom)

    # Other Functionality

    def RadiusOfGyration(self,chains=None,ismodel=False):
        
        ''' Acquire the radius of the gyration of the entire, or a portion of, the
        PDB protein molecule. '''
        
        if ismodel:
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
        
        ''' Populates the centroids of all residues. '''
        
        if not chain in self.chains:
            raise ValueError('Chain %s does not exist!' % (chain))

        out = []
        ch = self.chains[chain]

        for res in ch:
            res = ch[res]
            res.Centroid()
            out.append(res.centroid)

        return out
        
    def IndexSeq(self, chain, fst, ismodel=False):

        ''' Store in residues the correct index to the fasta.
        Requires a 1-to-1 correspondance at least a portion
        of the way through. Deprecated; use GetFASTAIndices(). '''
        
        if ismodel: thing = self.GetModel(chain)
        else: thing = self.GetChain(chain)
        return self.GetFASTAIndices(thing, fst)            
            
    def GetFASTAIndices(self, thing, fst):

        ''' Given a PDBchain, find 1-to-1 correspondances between
        it and a FASTA sequence. '''
        
        chainseq = thing.AsFASTA()
        ungapped = fst.replace('-','')
        if len(ungapped) != len(chainseq): return False
        for i in xrange(0,len(chainseq)):
            # See if there is a failed correspondance.
            if chainseq[i] != ungapped[i]: return False
        index = -1
        for i in thing:
            i = thing[i]
            index = fst.find(aa[i.name],index+1)
            i.fstindex = index
        return True        

    def IterAllResidues(self):
        
        ''' Produce an iterator to allow one to iterate over all possible residues. '''
        
        for ch in self.chains.keys():
            chain = self.GetChain(ch)
            for res in chain: yield chain[res]   
        for mo in self.models.keys():
            model = self.GetModel(mo)
            for res in model: yield model[res]

    def WriteGM(self,fasta,gm,chains=None,ismodel=False,CA=False):
        
        ''' Write the information present in this PDB between multiple
        chains as a Geometric Morphometric text file. This file will be
        formatted such that individual lines correspond to chains and semi-colons
        separate the (x,y,z) coordinates between all homologous residue positions. 
        Requires a FASTA alignment. Options include using alpha-carbon positions.
        By default, uses centroids of residues. '''
        
        if ismodel:
            orderofthings = self.orderofmodels
            things        = self.models
        else:
            orderofthings = self.orderofchains
            things        = self.chains
        if chains == None: chains = orderofthings
            
        fgm = open(gm,'w')
            
        # Acquire all homologous positions as defined in FASTA.
        mapng,masks,posvect = self.GetResidueAssociations(fasta,chains,ismodel)

        # Set up name grabbing from FASTA file for GM labelling.
        f = FASTAnet.FASTAstructure(fasta,uniqueOnly=False)
        name = lambda d: f.orderedSequences[d].name
        
        # Start writing.
        for ch in mapng:
            nm = name(mapng[ch]).split(':')[0]
            fgm.write('>%s:%s;' % (nm,nm[:4]))
            for pos in posvect:
                res = masks[mapng[ch]][pos]
                if CA: atom = res.GetCA()
                else:  atom = res.Centroid()
                fgm.write('%f;%f;%f;' % (atom.x,atom.y,atom.z))
            fgm.write('\n')
            
    def WriteLandmarks(self,fasta,lm,chains=None,ismodel=False):
        
        ''' Write the information present in this PDB between multiple
        chains as a landmark text file. This file will be formatted such that 
        the file is partitioned in sections starting with chain names and individual 
        lines in these correspond to homologous residue positions denoted by
        homologous position, residue number, and residue name tab-delimited. Requires
        a FASTA file. '''
        
        if ismodel:
            orderofthings = self.orderofmodels
            things        = self.models
        else:
            orderofthings = self.orderofchains
            things        = self.chains
        if chains == None: chains = orderofthings
            
        flm = open(lm,'w')
            
        # Acquire all homologous positions as defined in FASTA.
        mapng,masks,posvect = self.GetResidueAssociations(fasta,chains,ismodel)

        # Set up name grabbing from FASTA file for labelling.
        f = FASTAnet.FASTAstructure(fasta,uniqueOnly=False)
        name = lambda d: f.orderedSequences[d].name
        
        # Start writing.
        for ch in mapng:
            nm = name(mapng[ch]).split(':')[0]
            flm.write('>%s\n' % (nm))
            ind = 0
            for pos in posvect:
                res = masks[mapng[ch]][pos]
                flm.write('%s\t%s\t%s\n' % (ind,res.index,res.name))
                ind += 1

    def Contacts(self,chain='ALL',thres=4.5):
        
        ''' Compute the contact map of all chains or a chain. '''     
        
        done = []
        if chain != 'ALL':
            chain = self.GetChain(chain)
            for rA in chain:
                for rB in chain:
                    resA,resB = chain[rA],chain[rB]
                    if resA == resB or (resA,resB) in done: continue  
                    elif resA.InContactWith(resB,thres):
                        if not (int(resA),int(resB)) in self.contactmap:
                            self.contactmap.append((int(resA),int(resB)))
                            self.contactmap.append((int(resB),int(resA)))
                    done.append((resA,resB))
                    done.append((resB,resA))
        else: # All chains.
            for resA in self.IterAllResidues():
                for resB in self.IterAllResidues():
                    if resA == resB or (resA,resB) in done: continue 
                    elif resA.InContactWith(resB,thres):
                        A = resA.index+resA.chain
                        B = resB.index+resB.chain
                        if not (A,B) in self.contactmap:
                            self.contactmap.append((A,B))
                            self.contactmap.append((B,A))
                    done.append((resA,resB))
                    done.append((resB,resA))
                    
        self.contactmap = sorted(self.contactmap, key=lambda element: (
            element[0], element[1]))

    def WriteContacts(self, filename):
        
        ''' Write contact map. '''
        
        fout = open(filename, 'w')
        for a in self.contactmap:
            fout.write('%s\n'%(str(a)))
        fout.close()            

    # Internals

    def __iter__(self):
        
        ''' Returns all PDBchains as an iterator. '''
        
        for ch in self.chains: yield self.chains[ch]
        for mo in self.models: yield self.models[mo]

    def __len__(self):

        ''' Returns the length of the protein. '''

        # For all chains.
        chs = self.chains
        return sum([len(chs[x]) for x in chs])

    def _remarksToString(self):
        
        remarkstr = ''
        for it in xrange(len(self.remarks)):
            remarkstr += 'REMARK%4s %s\n' % (str(it+1),self.remarks[it])
        return remarkstr

    def __str__(self):

        ''' As a string, outputs structure in the PDB format. '''

        out = ''
        out += self._remarksToString()
        for chain in self.orderofchains: out += str(self.GetChain(chain))
        for model in self.orderofmodels: out += str(self.GetModel(model))
        return out

if __name__ == "__main__":
    mystructure = PDBstructure(sys.argv[1])
    if len(sys.argv) > 2:
        print 'RMSD',mystructure.rmsd(sys.argv[2])
        print 'RRMSD',mystructure.rrmsd(sys.argv[2])
        print 'TMscore',mystructure.tmscore(sys.argv[2])
        print 'GDT',mystructure.gdt(sys.argv[2])
    x =  mystructure.GetAllCentroid('A')
    mystructure.Contacts()
    print mystructure.contactmap
