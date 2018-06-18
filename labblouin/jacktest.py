from labblouin import PDBnet
import numpy as np
from MDgraph import MDgraph

s = MDgraph.ContactMatrixOperations(filein='/home/kestrel/locallab/pdb/2BEGmichelle.pdb', step=10,
                                    M=np.load('/home/kestrel/locallab/cm/2BEGmichellethres4.5step10-FreqCM.npy'))
# s.MakeFreqMatrix(step=10)
s.ComputeTriSuperMatrix(step=10)
# s.SaveFreqMatrix(savetarget='/home/kestrel/locallab/cm/',filename='2BEGmichelle', step=10)
s.ComputeTriSuperMatrix(step=10)
s.SaveSuperMatrix(savetarget='/home/kestrel/locallab/sm/', filename='2BEGmichelle', step=10)