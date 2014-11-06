from array import array
from collections import MutableMapping
from itertools import izip

# Constants

DICT_START_NUM  = 8
DICT_FREE_SLOT  = -1
DICT_DUMMY_SLOT = -2
DICT_SHIFT      = 5

class memEfficientDictionary(MutableMapping):
    
    ''' A memory-efficient alternative to the built-in Python dictionary with a trade-off
    for (insertion) performance. Based on recipe found by Raymond Hettinger
    (http://code.activestate.com/recipes/578375-proof-of-concept-for-a-more-space-efficient
    -faster/). '''
    
    def __init__(self,*args,**kwds):
        
        if not hasattr(self,'keylist'): self.clear()
        self.update(*args,**kwds)
        
    def __len__(self): return self.used
    def __iter__(self): return iter(self.keylist)
    def iterkeys(self): return __iter__()
    def keys(self): return list(self.keylist)
    def itervalues(self): return iter(self.valuelist)
    def values(self): return list(self.valuelist)
    def iteritems(self): return izip(self.keylist,self.valuelist)
    def items(self): return zip(self.keylist,self.valuelist)
    
    def __contains__(self,key):
        index,i = self._access(key,hash(key))
        return index >= 0
    
    def get(self,key,default=None):
        index,i = self._access(key,hash(key))
        if index >= 0: return self.valuelist[index]
        return default
    
    def popitem(self):
        assert self.keylist
        key = self.keylist[-1]
        val = self.valuelist[-1]
        del self[key]
        return key,val
        
    def clear(self):
        
        ''' Clear the data structure. '''
        
        self.indices   = self._index(DICT_START_NUM)
        self.hashlist  = []
        self.keylist   = []
        self.valuelist = []
        self.used      = 0
        self.filled    = 0 # Includes dummies.
    
    def __getitem__(self,key):
        
        ''' Access an item through the [] notation. '''
        
        hashv = hash(key)
        index,i = self._access(key,hashv)
        if index < 0: return None
        return self.valuelist[index]
    
    def __setitem__(self,key,value):
        
        hashv = hash(key)
        index,i = self._access(key,hashv)
        if index < 0: # Does not exist.
            self.indices[i] = self.used
            self.hashlist.append(hashv)
            self.keylist.append(key)
            self.valuelist.append(value)
            self.used += 1
            if index == DICT_FREE_SLOT:
                self.filled += 1
                if self.filled * 3 > len(self.indices) * 2:
                    self._resize(4 * len(self))
        else: self.valuelist[index] = value
    
    def __delitem__(self,key):
        
        hashv = hash(key)
        index,i = self._access(key,hashv)
        if index < 0: raise KeyError(key)
        self.indices[i] = DICT_DUMMY_SLOT
        self.used -= 1
        # Swap with last entry to avoid leaving hole.
        if index != self.used:
            lasth = self.hashlist[-1]
            lastk = self.hashlist[-1]
            lastv = self.valuelist[-1]
            lasti,j = self._access(lastk,lasth)
            self.indices[j] = index
            self.hashlist[index] = lasth
            self.keylist[index] = lastk
            self.valuelist[index] = lastv
        self.hashlist.pop()
        self.keylist.pop()
        self.valuelist.pop()
    
    @staticmethod
    def _index(n):
        
        if (n <= 2**7): return array('b',[DICT_FREE_SLOT])*n    # Signed Character
        elif (n <= 2**15): return array('h',[DICT_FREE_SLOT])*n # Signed Short
        elif (n <= 2**31): return array('l',[DICT_FREE_SLOT])*n # Signed Long
        return [DICT_FREE_SLOT] *n                              # Python Integers
    
    @staticmethod
    def _probe(hashv,mask):
        
        ''' Mimic current dictionary probing behavior. '''
        
        if hashv < 0: hashv = -hashv
        i = hashv & mask
        yield i 
        probe = hashv
        while 1:
            yield ((5*i+probe+1) & 0xFFFFFFFFFFFFFFFF) & mask
            probe >>= DICT_SHIFT
            
    def _access(self,key,hashv):
        
        ''' Lookup logic. '''
        
        freeSlot = None
        for i in self._probe(hashv,len(self.indices)-1):
            index = self.indices[i]
            if index == DICT_FREE_SLOT:
                if freeSlot is None: return (DICT_FREE_SLOT,i)
                else:                return (DICT_DUMMY_SLOT,freeSlot)
            elif index == DICT_DUMMY_SLOT:
                if freeSlot is None: freeSlot = i
            elif (self.keylist[index] is key or 
                  self.hashlist[index] == hashv and 
                  self.keylist[index] == key): return (index,i)
            
    def _resize(self,n):
        
        ''' Reindex all entries (hash, key, value) without moving
        entries. '''
        
        n = 2 ** n.bit_length()
        self.indices = self._index(n)
        for index,hashv in enumerate(self.hashlist):
            for i in self._probe(hashv,n-1):
                if self.indices[i] == DICT_FREE_SLOT: break
            self.indices[i] = index
        self.filled = self.used