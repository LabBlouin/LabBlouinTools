#!/usr/bin/python
'''
   Process MATT output and convert them to the GM format.
   
   MATT2GM Copyright (C) 2012 Christian Blouin
   MATT2GM Version Copyright (C) 2013 Christian Blouin - Jose Sergio Hleap - Alex Safatli
   
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
   
   E-mail: cblouin@cs.dal.ca, jshleap@squalus.org, iltafas@gmail.com
   
   Requirements:
   1) Python:
	  a) PDBnet.py
	  #######################################
	  #PDBnet.py is a python script and have#
	  #to be provided with this code        #
	  #######################################
   
'''

# Import section
import os
import sys
from PDBnet import PDBstructure
# Constants
aa = {'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE','G':'GLY','H':'HIS','I':'ILE','K':'LYS','L':'LEU','M':'MET','N':'ASN','P':'PRO','Q':'GLN','R':'ARG','S':'SER','T':'THR','V':'VAL','W':'TRP','Y':'TYR'}

# definition section
def LoadFasta(filename):
	# Load Fasta into a dictionary of sequences indexed by chain name
	data = open(filename).read().replace(':>',':####boo')
	data = data.split('>')[1:]
	for i in range(len(data)):
		if ':####boo' in data[i]:
			data[i] = data[i].replace(':####boo', ':>')

	# Output data structure
	out = {}
	outorder = []

	# Read in the sequences
	for i in data:
		# chain name
		colon = i.find(':')
		endline = i.find('\n')
		chain = i[colon+1:endline]

		# Sequence
		seq = i[endline:]
		seq = seq.replace('\n','')

		out[chain] = seq
		outorder.append(chain)

	return out, outorder

def GetMask(A):
	# Returns a string with * for all homologous sites and . for non all-homologous sites.
	# Assume that no gaps can occur in such sites

	# Output data structure
	out = ''

	# Get the lenght of any sequence (they are all the same length)
	L = len(A[A.keys()[0]])

	keys = A.keys()
	for i in range(L):
		c = '*'
		for k in keys:
			if A[k][i] == '-':
				c = '.'
				break
		out += c

	return out

def RemoveNonHomologous(pdb, aln, msk, mapfile):
	# Delete the residues that are not homologous through the whole alignment

	fout = open(mapfile,'w')

	# Get all chains
	keys = pdb.chains.keys()

	for chain in keys:
		struct = pdb.chains[chain]
		seq = aln[chain]
		indexlist = pdb.chains[chain].GetIndices()
		index = 0

		fout.write('>'+chain + '\n')
		counter = 0

		for i in indexlist:
			if msk[struct[i].fstindex] != '*':
				# Delete
				del struct[i]
			else:
				# Keep
				fout.write('%d\t%s\t%s\n'%(counter,struct[i].index.strip(),struct[i].name))
				counter += 1


	fout.close()

def Rename_fasta(prefix):
	'''
	Rename the chains in the fasta file
	'''
	# chain names
	chains = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z',
		      '1','2','3','4','5','6','7','8','9','0','~','!','@','#','$','%','^','&','(',')','-','+','=','_','<','>',
		      '/',"'","\\",'?','"','[',']','{','}','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q',
		      'r','s','t','u','v','w', 'x', 'y','z']
	# create a backup of the fasta
	os.system('cp %s.fasta ./%s.fasta_backup'%(prefix,prefix))
	#get the fasta file
	names=[]
	seqs=[]
	fasta=open(prefix+'.fasta').read().split('>')
	c=0
	for line in fasta[1:]:
		bline=line.split(':')
		names.append('>'+bline[0]+':'+chains[c])
		c+=1
		seqs.append(bline[1][bline[1].find('\n')+1:])

	#open outfile
	out=open(prefix+'.fasta','w')
	#write the new fastafile
	for j in range(len(names)):
		out.write(names[j]+'\n'+seqs[j])

	out.close()

def main(prefix,alphaC=False,redundant=False):
	# Form the proper file names
	pdbfile = prefix + '.pdb'
	fastafile = prefix + '.fasta'
	mapfile = prefix + '.landmarks'
	outfile = prefix + '.gm'

	# Read all structures into RAM
	pdbs = PDBstructure(pdbfile)

	# if redundant structure fasta need to rename the chains
	if redundant:
		Rename_fasta(prefix)

	# Get the alignment
	alns, kalns = LoadFasta(fastafile)

	# Check that both the alignment and the pdb have the same chain names
	kpdbs= pdbs.orderofchains
	newalns={}
	if not len(kalns) == len(kpdbs):
		print 'not the same number of chains in the PDB and FASTA files'
		sys.exit(-1)
	else:
		for a in range(len(kalns)):
			newalns[kpdbs[a]] = alns.pop(kalns[a])
	alns = None
	alns = newalns
				
	for chain in pdbs.chains:
		pdbs.IndexSeq(chain, alns[chain])

	# Get the mask
	mask = GetMask(alns)

	# Delete unaliagned sites
	RemoveNonHomologous(pdbs, alns, mask, mapfile)

	# Write to an output
	fout = open(outfile,'w')
	for chain in pdbs.chains:
		fout.write('>'+prefix+':'+chain+';')
		for index in pdbs.chains[chain].GetIndices():
			# deleted item
			if not pdbs.chains[chain].GetResidueByIndex(index):
				continue
			# Write XYZ
			pdbs.chains[chain][index].Centroid()
			if alphaC:
				ca= pdbs.chains[chain][index].atoms['CA']
			else:
				ca = pdbs.chains[chain][index].centroid

			for c in [ca.x,ca.y,ca.z]:
				fout.write('%f;'%(c))
		# End of line
		fout.write('\n')
		
# Application code
if __name__ == '__main__':
	# Read in the file prefix to process
	if len(sys.argv) < 2:
		print "Usage: MATT2GM.py prefix [options: -alphaC, -redundant]"
		sys.exit()
	alphaC=False
	redundant=False
	prefix= sys.argv[1]
	for arg in sys.argv[1:]:
		if arg == '-alphaC':
			alphaC=True
		elif arg == '-redundant':
			redundant=True
	
	main(prefix,alphaC,redundant)

