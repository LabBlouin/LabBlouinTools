#!/bin/python

''' 
A library to manage interfacing with the MDSA database raw files (http://dna.cs.byu.edu/mdsas/index.shtml).

mdsa Python Library / Oct 8, 2014 / Alex Safatli

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

E-mail: safatli@cs.dal.ca
Dependencies: IO, FASTAnet (LabBlouinTools)

'''

import os, FASTAnet, IO
from glob import glob as files

mdsaTypes = ['smart','balibase','oxbench','prefab']

class mdsaAlignment:

    ''' Model a single MDSA alignment. '''

    def __init__(self,finame):

        self.path      = finame
        self.name      = IO.getFileName(finame)
        self.fasta     = FASTAnet.FASTAstructure(self.path,uniqueOnly=False)

    def getObject(self): return self.fasta
    def getNames(self): return self.fasta.getSequenceNames()
    def getNumSequences(self): return len(self.getSequences())
    def getSequences(self): return self.fasta.getSequences()
    def getPath(self): return self.path
    
    def writeFASTA(self,fi,names=None):
        if not names:
            self.fasta.writeFile(fi)
            return fi
        fh = open(fi,'w')
        fh.write(self.getFASTAfor(names))
        fh.close()
        return fi
    
    def getFASTAfor(self,names):
        sequ = ''
        out = ''
        for name in names:
            sequ = ''
            seq = str(self.fasta.getSequenceByName(name))
            for i in xrange(0,len(seq),70): sequ += seq[i:i+70] + '\n'
            out += '>%s\n%s' % (name,sequ)
        return out
    
    def getSequenceLength(self): return self.getAlignmentLength()
    def getAlignmentLength(self):
        return len(self.fasta.sequences[self.fasta.sequences.keys()[0]])
    
class mdsaDatabase:
    
    def __init__(self,dbpath,traverse=True):
        self.path      = dbpath
        self.files     = {}
        if (traverse): self.traverse()
        
    def __iter__(self):
        for it in self.files:
            yield self.files[it]
            
    def traverse(self):
        if len(self.files) == 0:
            fis = files(os.path.join(self.path,'*'))
            for fi in fis:
                f = mdsaAlignment(fi)
                self.files[fi] = f
                    
    def getPath(self): return self.path
    def getFiles(self): return self.files
