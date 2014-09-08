#!/bin/python

''' 
An (incomplete) parsing library for output from the T-Coffee Expresso executable.

expresso Python Library / Summer 2013 / Alex Safatli

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
Dependencies: -

'''

class expressoParser:
    
    ''' Parses Expresso output. '''
    
    def __init__(self,fin):
        self.filein = fin
        self.lines = {}
        self.data = None
        self.__read__()
        
    def __read__(self):
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
        