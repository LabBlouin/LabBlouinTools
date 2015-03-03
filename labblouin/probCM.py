__author__ = 'Jack'

from labblouin import PDBnet
from matplotlib.pylab import *
import numpy


p = PDBnet.PDBstructure(filein='/home/kestrel/ryan/pdb_files/NPC1ch.pdb')


print 'protein 1 chain names: {}'.format(p.GetChainNames())


print 'protein 1 model names: {}'.format(p.GetModelNames())


print 'starting contacts'
print p.Contacts(1)
print 'end contacts'






