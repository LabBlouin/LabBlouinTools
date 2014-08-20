# homology.py
# -------------------------
# May 22, 2013; Alex Safatli
# -------------------------
# Homstrad Database API

import glob, os, IO
import re as R

def hasNoFolders(path):
    for fi in glob.glob(os.path.join(path,'*')):
        if os.path.isdir(fi): return False
    return True

class homstradFolder:
    def __init__(self,foldrname):
        self.path = foldrname
        self.name = IO.getFileName(foldrname)
        self.files = glob.glob(os.path.join(foldrname,'*'))
        self._pdb  = ''
        self.molecules = []
        self.sequences = {}
        self.__manifest__()
    def __manifest__(self):
        # Find PIR file.
        pirfile = None
        pdbfile = None
        for fi in self.files:
            if not pirfile and fi.endswith('.pir'): pirfile = fi
            if not pdbfile and fi.endswith('.pdb'): pdbfile = fi
        if not pirfile:
            raise IOError('%s did not possess PIR file.' % (self.path))
        if not pdbfile:
            raise IOError('%s did not possess PDB file.' % (self.path))
        self._pdb = pdbfile
        # Parse PIR file for sequences.
        pirhandle = open(pirfile)
        struct, seq = None, ''
        for line in pirhandle:
            if line.startswith('>'):
                if struct:
                    self.molecules.append(struct)
                    self.sequences[struct] = seq
                    struct = None
                    seq = ''
                struct = line.strip('>').split(';')[-1].strip('\n')
            elif line.startswith('structure'): pass
            elif line == '\n': pass
            else: seq += line.strip().strip('*')
        if struct:
            self.sequences[struct] = seq
            self.molecules.append(struct)
        pirhandle.close()
        # Parse PDB file for name sequence. Ensure consistency.
        if len(self.sequences) > 2:
            tmp = self.molecules
            self.molecules = []
            pdbhandle = open(pdbfile)
            s = ''
            line = pdbhandle.readline()
            while line.startswith('REMARK'):
                s += line.strip('\n')
                line = pdbhandle.readline()
            mols = R.findall('\S*\s*chain',s)
            for mol in mols: self.molecules.append(mol.split()[0]) 
            pdbhandle.close()
            if len(self.molecules) != len(self.sequences):
                if len(self.molecules) == 0: self.molecules = tmp
                else:
                    raise IOError('PDB file did not match PIR structure series for %s.' % (
                        self.path))        
    def getNames(self): return self.molecules
    def getSequences(self): return self.sequences
    def getFiles(self): return self.files
    def getPath(self): return self.path
    def getFASTA(self): return self.getFASTAfor(self.molecules)
    def writeFASTA(self,fi,names=None):
        if not names: names = self.sequences
        fh = open(fi,'w')
        fh.write(self.getFASTAfor(names))
        fh.close()
        return fi
    def getFASTAfor(self,names):
        sequ = ''
        out = ''
        for name in names:
            sequ = ''
            seq = self.sequences[name]
            for i in xrange(0,len(seq),70): sequ += seq[i:i+70] + '\n'
            out += '>%s\n%s' % (name,sequ)
        return out
    def getAlignedPDB(self): return self._pdb
    def getPDBs(self):
        for mol in self.molecules:
            yield os.path.join(self.path,'%s.atm' % (mol))
    def getPDBfor(self,name):
        if name in self.molecules:
            return os.path.join(self.path,'%s.atm' % (name))
        else: return None
    def getAlignmentLength(self):
        return len(self.sequences[self.sequences.keys()[0]])
    
    
class homstradDatabase:
    def __init__(self,dbpath,traverse=True):
        self.path      = dbpath
        self.folders   = {}
        self.failed    = []
        self.succeeded = []
        if (traverse): self.traverse()
    def traverse(self):
        if len(self.folders) == 0:
            folders = glob.glob(os.path.join(self.path,'*'))
            for folder in folders:
                if os.path.isdir(folder) and hasNoFolders(folder):
                    try:
                        f = homstradFolder(folder)
                        self.folders[folder] = f
                        self.succeeded.append(folder)
                    except IOError: self.failed.append(folder)
    def getPath(self): return self.path
    def getFolders(self): return self.folders
    def getFailedCount(self): return len(self.failed)
    def getSucceededCount(self): return len(self.succeeded)