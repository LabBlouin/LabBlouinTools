# scop.py
# -------------------------
# May 28, 2013; Alex Safatli
# -------------------------
# SCOP Database API

import glob, os, urllib, pickle, random

SCOP_VERSION = '1.75'
SCOP_PARSE_FILE = 'http://scop.mrc-lmb.cam.ac.uk'+\
    '/scop/parse/dir.%s.scop.txt_%s' % ('des',SCOP_VERSION)
SCOP_PDBDB_FILE = 'http://scop.mrc-lmb.cam.ac.uk'+\
    '/scop/parse/dir.%s.scop.txt_%s' % ('cla',SCOP_VERSION)
LOCAL_PARSE_FILE = 'directory.scop'
LOCAL_PDBDB_FILE = 'pdbdb.scop'
LOCAL_FLAT_FILE  = 'scop.pickl' 

class scopHierarchy:
    classes       = {}
    folds         = {}
    superfamilies = {}
    families      = {}
    domains       = {}
    species       = {}
    entries       = {}
    def __init__(self,cacheloc):
        self.cache = cacheloc
    def __downloadparse__(self):
        dest = os.path.join(self.cache,LOCAL_PARSE_FILE)
        try: dest, hdrs = urllib.urlretrieve(SCOP_PARSE_FILE,dest)
        except: raise IOError('Could not download <%r>.' % (SCOP_PARSE_FILE))
    def __downloadpdbdb__(self):
        dest = os.path.join(self.cache,LOCAL_PDBDB_FILE)
        try: dest, hdrs = urllib.urlretrieve(SCOP_PDBDB_FILE,dest)
        except: raise IOError('Could not download <%r>.' % (SCOP_PDBDB_FILE))        
    def __parsefile__(self,fi):
        fh = open(fi)
        for line in fh:
            if line.startswith('#'): continue
            l = line.strip().split('\t')
            if len(l) != 5: continue
            sunid, entry, sccs, name, desc = l
            it = scopItem(sunid,sccs,name,desc)
            if (entry == 'cl'):   self.classes[sunid] = it
            elif (entry == 'cf'): self.folds[sunid] = it
            elif (entry == 'sf'): self.superfamilies[sunid] = it
            elif (entry == 'fa'): self.families[sunid] = it
            elif (entry == 'dm'): self.domains[sunid] = it
            elif (entry == 'sp'): self.species[sunid] = it
            elif (entry == 'px'): self.entries[sunid] = it
        fh.close()
    def __pdbdbfile__(self,fi):
        fh = open(fi)
        for line in fh:
            if line.startswith('#'): continue
            l = line.strip().split('\t')
            if len(l) != 6: continue
            sid, pdb, dom, sccs, sunidpx, sunidli = l
            ids = sunidli.split(',')
            for i in ids:
                try:
                    if i.startswith('cl'): self.classes[i.strip('cl=')].addPDB(pdb)
                    elif i.startswith('cf'): self.folds[i.strip('cf=')].addPDB(pdb)
                    elif i.startswith('sf'): self.superfamilies[i.strip('sf=')].addPDB(pdb)
                    elif i.startswith('fa'): self.families[i.strip('fa=')].addPDB(pdb)
                    elif i.startswith('dm'): self.domains[i.strip('dm=')].addPDB(pdb)
                    elif i.startswith('sp'): self.species[i.strip('sp=')].addPDB(pdb)
                    elif i.startswith('px'): self.entries[i.strip('px=')].addPDB(pdb)
                except:
                    print 'sunid %s did not correspond to a parsed directory item.' % (i)
        fh.close()    
    def __readpickle__(self,fi):
        fh = open(fi)
        self.classes, self.folds, self.superfamilies, \
            self.families,self.domains, self.species, \
            self.entries = pickle.load(fh)
        fh.close()
    def __writepickle__(self,fi):
        fh = open(fi,'w')
        cl, fo, sf, ff, do, sp, ee = self.classes, self.folds, \
            self.superfamilies, self.families, \
            self.domains, self.species, self.entries
        pickle.dump((cl, fo, sf, ff, do, sp, ee),fh)
        fh.close()
    def populateHierarchy(self):
        fiflat = os.path.join(self.cache,LOCAL_FLAT_FILE)
        if not os.path.isfile(fiflat):
            fiparse = os.path.join(self.cache,LOCAL_PARSE_FILE)
            fipdbdb = os.path.join(self.cache,LOCAL_PDBDB_FILE)
            if not os.path.isfile(fiparse): self.__downloadparse__()
            if not os.path.isfile(fipdbdb): self.__downloadpdbdb__()
            self.__parsefile__(fiparse)
            self.__pdbdbfile__(fipdbdb)
            self.__writepickle__(fiflat)
        else: self.__readpickle__(fiflat)
    def getDissimilar(self,suprfquery,num=100):
        # Given a superfamily, grab a set of PDB files
        # from other superfamilies (1 random ambassador
        # from each).
        key = self.query(self.superfamilies,suprfquery)
        pdbs = []
        if not key: return None, None
        oth = [self.superfamilies[x] for x in self.superfamilies if x != key]
        random.shuffle(oth)
        oth = oth[:num]
        for it in oth:
            pdbli = [x for x in it.pdbs]
            random.shuffle(pdbli)
            pdbs.append(pdbli[0])
        # Get the PDB for the superfamily query.
        it = self.superfamilies[key]
        pdbli = [x for x in it.pdbs]
        random.shuffle(pdbli)
        sfpdb = pdbli[0]
        return sfpdb, pdbs
    def query(self,dicti,q):
        key = None
        if q in dicti: key = q
        else:
            for x in dicti:
                for it in [dicti[x].sccs, 
                           dicti[x].shortname, 
                           dicti[x].desc]:
                    if it.lower().find(q.lower()) != -1:
                        key = x
                        break
        return key       
    def getSimilar(self,suprfquery,num=100):
        # Given a superfamily, grab a set of PDB files
        # from that superfamilies.
        pdbs = []
        key = self.query(self.superfamilies,suprfquery)
        if not key: return None, None
        sf = self.superfamilies[key]
        pdbli = [x for x in sf.pdbs]
        random.shuffle(pdbli)
        pdbs = pdbli[0:num]
        if len(pdbs) > num+1:
            sfpdb = pdbli[num+1]
            return sfpdb, pdbs
        else: return None, pdbs
    def getFamilies(self,suprfquery):
        # Given a superfamily, grab a set of all associated
        # families.
        key = self.query(self.superfamilies,suprfquery)
        if not key: return []
        sf = self.superfamilies[key]
        sfpdbs = [x for x in sf.pdbs]
        families = []
        for fam in self.families:
            fpdbs = [x for x in self.families[fam].pdbs]
            if len(set(fpdbs).intersection(set(sfpdbs))) == len(fpdbs):
                families.append(self.families[fam])
        return families
    
class scopItem:
    def __init__(self,sunid,sccs,shortname,desc):
        self.sunid = sunid
        self.sccs = sccs
        self.shortname = shortname
        self.desc = desc
        self.pdbs = []
    def addPDB(self,pdb):
        if not pdb in self.pdbs: self.pdbs.append(pdb)
    def getPDBs(self):
        return self.pdbs

