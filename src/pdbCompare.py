#!/bin/python

# pdbCompare.py
# -------------------------
# May 16, 2012; Alex Safatli
# -------------------------
# Functions to help in comparing
# two sets of PDB files.
#
# Makes use of a special sort of md5
# signature in order to uniquely identify
# a PDB file by acquiring ATOM lines with
# alpha-carbons (CA).
#
# Can be run from command line with two
# arguments.

import pfam
import glob
import os
import md5
import sys
    
def getChecksum(pdb_file):
    '''
    Given a PDB file, retrieve a
    unique checksum.
    '''        
    if not os.path.isfile(pdb_file):
        return None
    f_in = open(pdb_file, 'r')
    checksum = md5.new()    
    for line in f_in:
        if 'ATOM' in line and ' CA ' in line:
            checksum.update(line.strip())
    return checksum.digest()
    
def pdbCrosscheck(folder1, folder2):
    '''
    Given two folder paths, determine what PDBs
    are unique amongst the two
    '''    
    if (not os.path.isdir(folder1)) or \
       (not os.path.isdir(folder2)):
        return None
    notunique = []   
    notunique_in_1 = []
    notunique_in_2 = []
    pdbList1 = glob.glob(os.path.join(folder1, '*.pdb'))
    pdbList2 = glob.glob(os.path.join(folder2, '*.pdb'))
    for pdb1 in pdbList1:
        for pdb2 in pdbList2:
            if getChecksum(pdb1) == getChecksum(pdb2):
                notunique.append(pdb1)
                notunique.append(pdb2)
    for pdb in pdbList1:
        if pdb in notunique:
            notunique_in_1.append(pdb)
    for pdb in pdbList2:
        if pdb in notunique:
            notunique_in_2.append(pdb)
    return (notunique_in_1,notunique_in_2, folder1, folder2)

def outputResults(pdbcrosscheck_results, fpath=None):
    compare = pdbcrosscheck_results
    str_out = 'All non-unique PDBs in %s:\n%s\n' % (compare[2], \
                                            pfam.list2txttable(compare[0], compare[2]))
    str_out += 'All non-unique PDBs in %s:\n%s\n' % (compare[3], \
                                            pfam.list2txttable(compare[1], compare[3]))
    if fpath is not None:
        fout = open(fpath, 'w')
        fout.write(str_out)
        fout.close()
    return str_out

# ----------------- If run from command-line -------------------- #

if __name__ == "__main__":

    if len(sys.argv) != 3:
        print 'Usage: pdbCompare.py folder1 folder2\n'
        exit()
    
    folder1 = sys.argv[1]
    folder2 = sys.argv[2]
    
    compare = pdbCrosscheck(folder1, folder2)
    
    if compare is None:
        raise SystemError('Error encountered with crosscheck.')
        exit()
    
    for c in xrange(len(compare[0])):
        compare[0][c] = os.path.splitext(os.path.basename(compare[0][c]))[0]
    for c in xrange(len(compare[1])):
        compare[1][c] = os.path.splitext(os.path.basename(compare[1][c]))[0]
    
    print outputResults(compare)
