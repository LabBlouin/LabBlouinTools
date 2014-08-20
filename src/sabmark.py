# sabmark.py
# -------------------------
# Jun 14, 2013; Alex Safatli
# -------------------------
# SABMARK Database API

import os
from glob import glob as g
from utils import IO

class sabmarkFile:
    def __init__(self,path,l,fp):
        self.path = path
        self.length = l
        self.true = fp

class sabmarkSummary:
    def __init__(self,path):
        self.path = path
        self.entries = []
        self.names = []
        self.__read__()
    def __read__(self):
        f = open(self.path)
        l = f.read()
        f.close()
        for x in l.split("\n")[1:]: 
            s = x.split()
            if len(s) <= 0: continue
            self.entries.append(s)
            self.names.append(s[0])

class sabmarkFolder:
    def __init__(self,foldrname):
        self.path = foldrname
        self.files = g(os.path.join(foldrname,'*'))
        self.name = ""
        self.pdbs = []
        self.pdbs_nf = []
        self.__manifest__()
    def __manifest__(self):
        grsums = os.path.join(self.path,"group.summary")
        self.name = IO.getFileName(self.path)
        pdbs = os.path.join(self.path,"pdb")
        if (not os.path.isfile(grsums)): raise IOError("Group summary could not be found.")
        elif (not os.path.isdir(pdbs)):  raise IOError("PDB directory for group not found.")
        summary = sabmarkSummary(grsums)
        for name in summary.names:
            e = summary.names.index(name)
            entry = summary.entries[e]
            length = int(entry[1].strip())
            fp = entry[2].strip() == '1' # is true?
            fpath = os.path.join(pdbs,name+'.ent')
            if not os.path.isfile(fpath):
                self.pdbs_nf.append(sabmarkFile(fpath,length,fp))
                print 'Warning: %s could not be found but existed in summary.' % (name)
            else: self.pdbs.append(sabmarkFile(fpath,length,fp))
        
class sabmarkDatabase:
    def __init__(self,dbpath):
        self.path = dbpath
        self.categories = ['sup','sup_fp','twi','twi_fp']
        self.folders = {}
        self.pdbs = []
        self.falsepositives = []
        self.truepositives  = []
    def traverse(self):
        for cat in [x for x in self.categories]:
            path = os.path.join(self.path,cat)
            if (not os.path.isdir(path)): self.categories.remove(cat)
            else: self.folders[(path,cat)] = []
        for path, cat in self.folders:
            groups = g(os.path.join(path,'group*'))
            if ('fp' in cat): self.falsepositives.append((path,cat))
            else: self.truepositives.append((path,cat))
            for group in groups:
                groupf = sabmarkFolder(group)
                self.folders[(path,cat)].append(groupf)
                self.pdbs.extend(groupf.pdbs)
    def groups(self,fkey):
        if (not fkey in self.folders): yield ('','')
        else:
            for group in self.folders[fkey]: yield group