#!/bin/python

''' An (incomplete) parsing library for output from the T-Coffee Expresso executable.
expresso Python Library / Summer 2013 / Alex Safatli

E-mail: safatli@cs.dal.ca
Dependencies: - '''

class expressoParser:
    
    ''' Parses Expresso output. '''
    
    def __init__(self,fin):
        self.filein = fin
        self.lines = {}
        self.data = None
        self.__read__()
        
    def __read__(self):
        
        ''' Read an Expresso output file; store the data inside this object. '''
        
        fh = open(self.filein)
        self.data = fh.readlines()
        fh.close()
        for li in self.data:
            splitted = li.split()
            if len(splitted) > 1: name, _, dat = splitted
            elif len(splitted) == 1: 
                name = splitted[0]
                dat = ''
            else: continue
            self.lines[name.strip('>')] = dat.strip()
            
    def getPDBcodes(self):
        
        ''' Return PDB code and chain as a tuple for each line. '''
        
        return [(x[:-1],x[-1]) for x in self.lines.values() if x]
    
    def writeNonEmpty(self,fout):
        
        ''' Writes non-empty PDB IDs to file. '''
        
        fh = open(fout,'w')
        for name in self.getNonEmptySequenceNames():
            fh.write('>%s _P_ %s\n' % (name,self.lines[name]))
        fh.close()
        
    def getNonEmptySequenceNames(self):
        
        ''' Gets non-empty PDB names. '''
        
        return [x for x in self.lines if x]
        