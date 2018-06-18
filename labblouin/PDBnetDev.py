#!/bin/python

"""
PDBnet is a collection of Python objects intended to model and contain
PDB protein data.

PDBnet Copyright (C) 2012 Christian Blouin
Contributions by Alex Safatli, Jose Sergio Hleap, and Jack Ryan
Major Refactoring (2014) done by Alex Safatli and Jose Sergio Hleap

E-mail: cblouin@cs.dal.ca, safatli@cs.dal.ca, jshleap@dal.ca
Dependencies:
	Scipy
	BioPython
	FASTAnet (contained in LabBlouinTools)
	intervaltree (pip install intervaltree)
	joblib (pip install joblib)
"""

from mmap import mmap as memoryMappedFile
import numpy as np
import sys, FASTAnet, math
from math import sqrt, acos, degrees
from os.path import getsize as sizeOfFile
from Bio.SeqUtils.ProtParam import ProteinAnalysis as PA
from scipy.spatial.distance import cdist
from joblib import Parallel, delayed
import multiprocessing
import warnings
from copy import deepcopy
import sklearn.decomposition as skl
import tempfile
from joblib import load, dump
import operator

import pyximport;

pyximport.install()
import distance  # This will show up as error in IDE, but it is a cython module
import itertools

# Metadata

__author__ = ['Christian Blouin', 'Jose Sergio Hleap', 'Alexander Safatli', 
              'Jack Ryan', 'Michelle Lu']
__version__ = '1.1.1'

# Constants

PDB_LARGE_FILE_SIZE = 100000000  # bytes

# Classes

# Auxiliary methods
def numonly(index):
    ''' Given a string, return only the numerical part'''
    newi = ''
    for i in index:
        if i.isdigit(): newi += i
    return newi


def singlerun(state, thres, numstates):

    '''
    PDBstructure External Method
    Multiprocessing for AggregateMatrix() method
    '''

    l = len(state)

    m = np.zeros((l, l))
    contact_tuples = state.ContactMap(thres=thres)
    print 'Progress: {}/{} states'.format(state.name, numstates)
    for pair in contact_tuples:
        if pair[0] < l and pair[1] < l:
            # add the found contacts into the matrix
            m[pair[0]][pair[1]] = 1
            m[pair[1]][pair[0]] = 1

    # return m
    return m, state.name


# compute Theobald's quaterion-based characteristic polynomial
# taken from https://gist.github.com/rmcgibbo

def _rmsd_qcp(conformation1, conformation2):
    """
    PRIVATE!!
    Compute the RMSD with Theobald's quaterion-based characteristic
    polynomial

    Rapid calculation of RMSDs using a quaternion-based characteristic polynomial.
    Acta Crystallogr A 61(4):478-480.

    THIS CODE IS TAKEN FROM https://gist.github.com/rmcgibbo

    Parameters
    ----------
    conformation1 : np.ndarray, shape=(n_atoms, 3)
    The cartesian coordinates of the first conformation
    conformation2 : np.ndarray, shape=(n_atoms, 3)
    The cartesian coordinates of the second conformation
    Returns
    -------
    rmsd : float
    The root-mean square deviation after alignment between the two pointsets
    """

    def _center(conformation):
        """
        Center and typecheck the conformation
        """

        conformation = np.asarray(conformation)
        if not conformation.ndim == 2:
            raise ValueError('conformation must be two dimensional')
        _, three = conformation.shape
        if not three == 3:
            raise ValueError('conformation second dimension must be 3')

        centroid = np.mean(conformation, axis=0)
        centered = conformation - centroid
        return centered

    # center and typecheck the conformations
    A = _center(conformation1)
    B = _center(conformation2)
    if not A.shape[0] == B.shape[0]:
        raise ValueError('conformation1 and conformation2 must have same number of atoms')
    n_atoms = len(A)

    # the inner product of the structures A and B
    G_A = np.einsum('ij,ij', A, A)
    G_B = np.einsum('ij,ij', B, B)

    # M is the inner product of the matrices A and B
    M = np.dot(B.T, A)

    # unpack the elements
    Sxx, Sxy, Sxz = M[0, :]
    Syx, Syy, Syz = M[1, :]
    Szx, Szy, Szz = M[2, :]

    # do some intermediate computations to assemble the characteristic
    # polynomial
    Sxx2 = Sxx * Sxx
    Syy2 = Syy * Syy
    Szz2 = Szz * Szz

    Sxy2 = Sxy * Sxy
    Syz2 = Syz * Syz
    Sxz2 = Sxz * Sxz

    Syx2 = Syx * Syx
    Szy2 = Szy * Szy
    Szx2 = Szx * Szx

    SyzSzymSyySzz2 = 2.0 * (Syz * Szy - Syy * Szz)
    Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2

    # two of the coefficients
    C2 = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2)
    C1 = 8.0 * (
        Sxx * Syz * Szy + Syy * Szx * Sxz + Szz * Sxy * Syx - Sxx * Syy * Szz - Syz * Szx * Sxy - Szy * Syx * Sxz)

    SxzpSzx = Sxz + Szx
    SyzpSzy = Syz + Szy
    SxypSyx = Sxy + Syx
    SyzmSzy = Syz - Szy
    SxzmSzx = Sxz - Szx
    SxymSyx = Sxy - Syx
    SxxpSyy = Sxx + Syy
    SxxmSyy = Sxx - Syy
    Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2

    # the other coefficient
    C0 = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2 \
        + (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2) \
        + (-(SxzpSzx) * (SyzmSzy) + (SxymSyx) * (SxxmSyy - Szz)) * (
            -(SxzmSzx) * (SyzpSzy) + (SxymSyx) * (SxxmSyy + Szz)) \
        + (-(SxzpSzx) * (SyzpSzy) - (SxypSyx) * (SxxpSyy - Szz)) * (
            -(SxzmSzx) * (SyzmSzy) - (SxypSyx) * (SxxpSyy + Szz)) \
        + (+(SxypSyx) * (SyzpSzy) + (SxzpSzx) * (SxxmSyy + Szz)) * (
            -(SxymSyx) * (SyzmSzy) + (SxzpSzx) * (SxxpSyy + Szz)) \
        + (+(SxypSyx) * (SyzmSzy) + (SxzmSzx) * (SxxmSyy - Szz)) * (
            -(SxymSyx) * (SyzpSzy) + (SxzmSzx) * (SxxpSyy - Szz))

    # Netwtorn-Raphson
    E0 = (G_A + G_B) / 2.0
    max_eigenvalue = E0
    for i in range(50):
        old_g = max_eigenvalue
        x2 = max_eigenvalue * max_eigenvalue
        b = (x2 + C2) * max_eigenvalue
        a = b + C1
        delta = ((a * max_eigenvalue + C0) / (2.0 * x2 * max_eigenvalue + b + a))
        max_eigenvalue -= delta
        if abs(max_eigenvalue - old_g) < abs(1e-11 * max_eigenvalue):
            break
    if i >= 50:
        raise ValueError('More than 50 iterations needed.')

    rmsd = np.sqrt(np.abs(2.0 * (E0 - max_eigenvalue) / n_atoms))
    return rmsd


class PDBatom(object):
    """ ATOM in a PDB protein structure. """

    __slots__ = ('name', 'serial', 'x', 'y', 'z', 'occupancy', 'tempFactor',
                 'charge', 'symbol', 'parent')

    def __init__(self, serial, name, x, y, z, oc, b, symbol, charge):

        """ Construct an atom in a PDB protein structure. """

        self.name = name
        self.serial = serial
        self.x = x
        self.y = y
        self.z = z
        self.occupancy = oc
        self.tempFactor = b
        self.charge = charge
        self.symbol = symbol
        self.parent = None

    def fixname(self):

        """ Ensures the name of this atom fits within 4 characters. """

        if len(self.name) == 4:
            return self.name
        else:
            return ' ' + self.name.ljust(3)

    def __str__(self):

        """ 
        Express this atom as a single line in a PDB file (see PDB format). 
        """
        line = "ATOM  %5s %4s %3s %1s%4s    %8.3f%8.3f%8.3f%6.2f%6.2f       "\
            "   %2s%2s"
        return  line % (self.serial, self.fixname(), self.parent.name, 
                        self.parent.chain, self.parent.index, self.x, self.y, 
                        self.z, self.occupancy, self.tempFactor, self.symbol, 
                        self.charge)

    def Element(self):

        """ Get the atom names. """

        return (self.fixname().strip())

    def ThreeLetters(self):

        """ Get the residue names. """

        return (self.parent.name)

    def GetPosition(self):

        """ Get the 3-dimensional coordinate of this atom in space. """

        return (self.x, self.y, self.z)

    def DistanceTo(self, atom):

        """ Acquire the distance from this atom to another atom. """

        # return ((self.x-atom.x)**2+(self.y-atom.y)**2+(self.z-atom.z)**2)**0.5

        return distance.dist(self.GetPosition(), atom.GetPosition())

    def Bond(self, atom):

        """ Get the bond vector between two atoms.

		:param atom: the atom that is on the other side of this bond
		"""

        vector = (self.x - atom.x, self.y - atom.y, self.z - atom.z)

        return vector


class PDBterminator(PDBatom):
    """ 
    A placeholder class that represents a terminating ATOM-like line in the 
    PDB file. 
    """

    __slots__ = ('name', 'serial', 'x', 'y', 'z', 'occupancy', 'tempFactor',
                 'charge', 'symbol', 'parent', 'lastatom', 'lastresname',
                 'lastreschain', 'lastresind')

    def __init__(self, chaininst):
        self.lastatom = chaininst.GetAtoms()[-1]
        self.lastresname = self.lastatom.parent.name
        self.lastreschain = self.lastatom.parent.chain
        self.lastresind = self.lastatom.parent.index
        newserial = str(int(self.lastatom.serial) + 1)
        super(PDBterminator, self).__init__(newserial, '', 0, 0, 0, 0, 0, '', 
                                            '')

    def __str__(self):
        return 'TER   %5s %4s %3s %1s%4s' % (
            self.serial, '', self.lastresname, self.lastreschain, 
            self.lastresind)


class PDBresidue:
    """ A residue (collection of ATOM fields) in a PDB protein structure. """

    __slots__ = ('index', 'name', 'atoms', 'atomsOrdered', 'contacts', 'centroid',
                 'chain', 'model', 'fstindex')

    def __init__(self, index=None, name='', makeboundingbox=False):

        """ Construct a residue (collection of atoms) in a PDB protein structure. """

        self.index = index
        self.name = name
        self.atoms = {}
        self.atomsOrdered = []
        self.contacts = []
        self.centroid = None
        self.chain = ''
        self.model = None
        # self.fstindex = None
        self.boundingbox = None
        if makeboundingbox: self._ComputeBoudingBox()

    def GetAtoms(self):
        return [x for x in self.atomsOrdered]

    def AddAtom(self, atom):

        """ Add a PDBatom structure to this residue. """

        self.atoms[atom.name.strip()] = atom
        self.atomsOrdered.append(atom)
        atom.parent = self

    def GetCA(self):

        """ Get the alpha-carbon found in this residue as a PDBatom. """

        if 'CA' in self.atoms:
            return self.atoms['CA']
        else:
            raise LookupError("Residue %s (%s, chain: %s) did not possess an alpha-carbon." % (
                self.index, self.name, self.chain))

    def Centroid(self):

        """ Calculate the centroid of this residue. Return this as a PDBatom. """

        x = 0.0
        y = 0.0
        z = 0.0

        for atom in self.atoms:
            if atom in ['C', 'N', 'O']: continue
            x += self.atoms[atom].x
            y += self.atoms[atom].y
            z += self.atoms[atom].z

        div = float(len(self.atoms) - 3)

        x /= div
        y /= div
        z /= div

        self.centroid = PDBatom(None, 'centroid', x, y, z, 0, 0, '', '')
        self.centroid.parent = self

        return self.centroid

    def BondAngle(self, vector1, vector2):

        """ Calculate a bond angle between two given bond vectors

		:param vector1: you can generate a bond vector from PDBatom.Bond
		:param vector2: same as above
		"""

        dot_product = vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2]

        length1 = distance.dist(vector1, (0.0, 0.0, 0.0))
        length2 = distance.dist(vector2, (0.0, 0.0, 0.0))

        angle = 180 - (math.acos(dot_product / (length1 * length2)) * 180/ math.pi)

        return angle

    def HydrogenBond(self, other, distance_cutoff=3.5, angle_cutoff=30):
        '''
        Detect hydrogen bonds (HBs) between given residues, default thresholds for bond length and bond angle of a HB are determined based on
        the paper "Geometrical Preferences of the Hydrogen Bonds of Protein-Ligand Binding Interface Derived from Statistical
        Surveys and Quantum Mechanics Calculations" (2008), please cite the paper if you need to include these information.
        A hydrogen bond will be detected if the distance between the donor and acceptor is less than 3.5 angstroms and the angle is larger
        than 90 degrees.

        :param other: another residue that you would like to detect if there is a hydrogen bond formed with
        :param distance_cutoff: for now the distance threshold has been set to 3.5 angstroms
        :param angle_cutoff: angle threshold has been set to 90 degrees
        '''

        # backbone hydrogen bonds detection

        hydrogen_bonds = {}

        for atom1 in self.atoms:
            if atom1 in ['H']:
                for atom2 in other.atoms:
                    if atom2 in ['O']:
                        d = self.atoms[atom1].DistanceTo(other.atoms[atom2])
                        vector1 = self.atoms[atom1].Bond(other.atoms[atom2])
                        vector2 = self.atoms[atom1].Bond(self.atoms['N'])
                        angle = self.BondAngle(vector1, vector2)
                        if d <= distance_cutoff and angle <= angle_cutoff:
                            hydrogen_bonds[(self.atoms[atom1].ThreeLetters() + '-main',
                                            other.atoms[atom2].ThreeLetters() + '-main')] = d


        # side chain hydrogen bonds detecton
        '''for atoms1 in self.atoms:
            if atoms1 in ['HD1', 'HD2', '1HD2', '2HD2', 'HE', 'HE1', 'HE2', '1HE2', '2HE2' 'HZ', 'HZ1', 'HZ2',
                          'HZ3' 'NH1',
                          '1HH1', '2HH1', '1HH2', '2HH2', 'OH', 'HG', 'HH']:
                for atoms2 in other.atoms:
                    if atoms2 in ['ND1', 'OD1', 'OD2', 'OE1', 'OE2', 'OG', 'SD']:
                        d = self.atoms[atom1].DistanceTo(other.atoms[atom2])
                        vector1 = self.atoms[atom1].Bond(other.atoms[atom2])
                        for a in self.atoms:
                            if self.atoms[str(a)].serial == self.atoms[atom1].serial - 1:
                                vector2 = self.atoms[atom1].Bond(self.atoms[a])
                        angle = self.BondAngle(vector1, vector2)
                        if d <= distance_cutoff and angle <= angle_cutoff:
                            hydrogen_bonds[(self.atoms[atom1].ThreeLetters() + '-side',
                                            other.atoms[atom2].ThreeLetters() + '-side')] = d'''
        if len(hydrogen_bonds) == 0:
            return None
        else:
            return hydrogen_bonds

    def InContactWith(self, other, thres=4.5):

        """ Determine if in contact with another residue. """

        d = self.GetCA().DistanceTo(other.GetCA())

        if d <= thres:  # The CA's are in contact
            return True
        elif d < 12:  # They COULD be in contact...
            for atomA in self.atoms:
                if atomA in ['C', 'N', 'O']:
                    continue
                atomA = self.atoms[atomA]

                for atomB in other.atoms:
                    if atomB in ['C', 'N', 'O']:
                        continue
                    atomB = other.atoms[atomB]

                    if atomA.DistanceTo(atomB) <= thres:
                        return True
        else:  # They're not in contact
            return False

    def _ComputeBoundingBox(self):
        '''
        Compute the bounding box in x, y, and z for this residue.
        '''

        self.boundingbox = {}

        vals = [atom.GetPosition() for atom in self.GetAtoms()]
        xvals = [x[0] for x in vals]
        yvals = [y[1] for y in vals]
        zvals = [z[2] for z in vals]

        min_x = min(xvals)
        max_x = max(xvals)
        min_y = min(yvals)
        max_y = max(yvals)
        min_z = min(zvals)
        max_z = max(zvals)

        self.boundingbox['minx'] = min_x
        self.boundingbox['miny'] = min_y
        self.boundingbox['minz'] = min_z

        self.boundingbox['maxx'] = max_x
        self.boundingbox['maxy'] = max_y
        self.boundingbox['maxz'] = max_z

    def __int__(self):
        return int(self.index)

    def __str__(self):
        return '\n'.join([str(x) for x in self.atomsOrdered])


class PDBchain(object):
    """ A PDB chain (collection of protein residues). """

    __slots__ = ('name', 'structure', 'residues', 'indices', 'parent', 'fstdict')

    def __init__(self, name):

        self.name = name
        self.residues = list()
        self.indices = list()
        self.structure = None
        self.parent = None
        self.fstdict = {}

    def GetResidues(self):

        """ Acquire all of the residues in this chain. Note that if this
		is a large file, this will force the object to load all of these
		residues into memory from the file. """

        return [self.GetResidueByIndex(i) for i in self.GetIndices()]

    def IterResidues(self):

        """ Iteratively yield all of the residues in this chain. Note that if
		this is a large file, this will force the object to load all of these
		residues into memory from the file. """

        for i in self.GetIndices(): yield self.GetResidueByIndex(i)

    def GetResidueByIndex(self, index):

        """
		Alias for acquiring a residue from this instance using [] operator.
		"""

        return self.__getitem__(index)

    def GetResidueByPosition(self, pos):

        return self.residues[pos]

    def GetAtoms(self):

        """ Get all comprising atoms in this chain in order of residues. """

        atoms = []
        for x in self.IterResidues():
            atoms.extend(x.GetAtoms())
        return atoms

    def GetIndices(self):

        """ Acquire all of the indices for residues present in this chain. """

        return [int(x) if type(x) == int else x for x in self.indices]

    def AddIndexOfResidue(self, index):

        """ Add an index for a residue to this chain. Will generate own residue
		information if called upon (forces lazy evaluation). """

        self.indices.append(str(index))

    def AddResidue(self, resid):

        """ Add a residue to this chain. """

        resid.chain = self.name
        self.residues.append(resid)
        self.indices.append(str(resid.index))

    def AddResidueByIndex(self, index):

        """ Add a residue to this chain by its index in the PDB. The residue
		object will automatically be constructed from the file. """

        self.indices.append(str(index))
        resid = self.GetResidueByIndex(index)
        self.residues.append(resid)

    def RemoveResidue(self, resid):

        """ Remove a residue from this chain. """

        return self.pop(resid)

    def PCAContactMap(self):

        contactmap = []

        reslist = self.GetResidues()

        vals = np.zeros((self.__len__(), 6))

        rc = [[res.boundingbox['minx'], res.boundingbox['maxx'],
               res.boundingbox['miny'], res.boundingbox['maxy'],
               res.boundingbox['minz'], res.boundingbox['maxz']] for res in self.GetResidues()]

        for i in xrange(len(rc)):
            vals[i][0] = rc[i][0]
            vals[i][1] = rc[i][1]
            vals[i][2] = rc[i][2]
            vals[i][3] = rc[i][3]
            vals[i][4] = rc[i][4]
            vals[i][5] = rc[i][5]

        rc = np.array(rc)

        pca = skl.PCA(1)

        fit = pca.fit_transform(vals)

        print pca.explained_variance_ratio_

        combolist = [c for c in itertools.combinations(range(self.__len__()), 2)]

        for combo in combolist:
            if abs(fit[combo[0]] - fit[combo[1]]) < 10:  # As of July 27, this value changes with each protein :(
                if reslist[combo[0]].InContactWith(reslist[combo[1]]):
                    contactmap.append(combo)

        return contactmap

    def BBContactMap(self, thres=4.5):

        '''
		Compute the contactmap of this chain, using bounding boxes and an intervaltree.
		'''

        contactmap = []

        if len(self.indices) != 0:  # Check for if a model (with a chain) has been passed in

            reslist = self.GetResidues()

            # Make three interval trees: one for each of x, y, and z coordinates for the min/max bounding boxes
            tx = IntervalTree(
                [Interval(res.boundingbox['minx'], res.boundingbox['maxx'], res.index) for res in reslist])
            ty = IntervalTree(
                [Interval(res.boundingbox['miny'], res.boundingbox['maxy'], res.index) for res in reslist])
            tz = IntervalTree(
                [Interval(res.boundingbox['minz'], res.boundingbox['maxz'], res.index) for res in reslist])

            for residue in reslist:
                # Remove the current residue from each of the trees

                resx = Interval(residue.boundingbox['minx'], residue.boundingbox['maxx'], residue.index)
                resy = Interval(residue.boundingbox['miny'], residue.boundingbox['maxy'], residue.index)
                resz = Interval(residue.boundingbox['minz'], residue.boundingbox['maxz'], residue.index)

                tx.remove(resx)
                ty.remove(resy)
                tz.remove(resz)

                slicex = set([x.data for x in sorted(tx[resx.begin - 4.5: resx.end + 4.5])])
                slicey = set([x.data for x in sorted(ty[resy.begin - 4.5: resy.end + 4.5])])
                slicez = set([x.data for x in sorted(tz[resz.begin - 4.5: resz.end + 4.5])])

                # Find the intersection of the three slices in x/y/z
                finalreslist = set.intersection(slicex, slicey, slicez)

                # Check if the filtered residues are actually in contact with the residue in question
                finalcontactslist = [(residue.index, res)
                                     for res in finalreslist if residue.InContactWith(self.GetResidueByIndex(res))]

                # Finally, extend the contactmap list to contain the index of the residues
                contactmap.extend(finalcontactslist)

        contactmap = [contact for contact in contactmap if contact[0] != contact[1]]

        return contactmap

    def CAContactMap(self):

        '''
		CAContactMap uses a coarse alpha-carbon comparison to determine if residue combinations should
		be investigated further for contact identity.
		'''

        contactmap = []

        if self.residues != 0:

            reslist = self.residues

            dlist = [(reslist.index(c[0]), reslist.index(c[1]),
                      distance.dist([c[0].GetCA().x, c[0].GetCA().y, c[0].GetCA().z],
                                    [c[1].GetCA().x, c[1].GetCA().y, c[1].GetCA().z]))
                     for c in itertools.combinations(reslist, 2)]

            darray = np.recarray(len(dlist), dtype=[('x', int), ('y', int), ('distance', float)])

            for i in range(len(dlist)):
                darray[i]['x'] = dlist[i][0]
                darray[i]['y'] = dlist[i][1]
                darray[i]['distance'] = dlist[i][2]

            # maybearray -- holds values between 4.5 and 9 (they MAY be in contact)
            maybearray = darray[np.logical_and(darray['distance'] > 4.5, darray['distance'] <= 9)]
            # surearray -- holds values below 4.5 (they ARE in contact_
            surearray = darray[darray['distance'] <= 4.5]

            for t in surearray:
                contactmap.append((t['x'], t['y']))

            for t in maybearray:
                resA = reslist[t['x']]
                resB = reslist[t['y']]
                if resA.InContactWith(resB):
                    contactmap.append((t['x'], t['y']))

        return contactmap

    def ContactMap(self, thres=4.5, userindices=None):

        """ Compute the contact map of this chain. """

        if userindices:
            Range = userindices
        else:
            Range = range(len(self.GetIndices()))

        contactmap = []

        '''
		# for rA in self:
		for rA in residues:
			# for rB in self:
			for rB in residues:
				# ADDED BY JACK @ rev 1181
				if rA == rB or (rA,rB) in done: continue
				elif rA.InContactWith(rB, thres):
					if not (int(rA),int(rB)) in contactmap:
						contactmap.append((rA.index,rB.index))
						contactmap.append((rB.index,rA.index))
				done[(rA,rB)] = True
				done[(rB,rA)] = True
				# TAKEN OUT BY JACK @ rev 1181

				resA,resB = self[rA],self[rB]
				if resA == resB or (resA,resB) in done: continue
				elif resA.InContactWith(resB,thres):
					if not (int(resA),int(resB)) in contactmap:
						contactmap.append((int(resA),int(resB)))
						contactmap.append((int(resB),int(resA)))
				done[(resA,resB)] = True
				done[(resB,resA)] = True

			return sorted(contactmap, key=lambda e: (e[0], e[1]))
		'''

        # Make all possible combinations of residue indices as tuples
        combos = itertools.combinations(Range, 2)
        reslist = self.GetResidues()

        for combo in combos:
            rA = reslist[combo[0]]
            rB = reslist[combo[1]]
            if rA.InContactWith(rB, thres):
                contactmap.append(combo)

        return contactmap

    def GetPrimaryPropertiesFromBioPython(self):

        """ Use BioPython to populate this class with attributes. """

        props = PA(self.AsFASTA())
        self.amino_acid_composition = props.get_amino_acids_percent()
        self.molecular_weight = props.molecular_weight()
        self.aromaticity = props.aromaticity()
        self.instability_index = props.instability_index()
        self.flexibility = props.flexibility()
        self.isoelectric_point = props.isoelectric_point()
        self.secondary_structure_fraction = props.secondary_structure_fraction()
        return props

    def AsFASTA(self):

        """ Return the string representing this chain's FASTA sequence. """

        return ''.join([aa[item.name] for item in self.IterResidues()])

    def WriteAsPDB(self, filename):

        """ Write this single chain to a file as a PDB. """

        fh = open(filename, 'w')
        fh.write(str(self))
        fh.close()

    def SortByNumericalIndices(self):

        """ Sort all internal items by a numerical index. """

        if self.structure.handle:
            if not self.structure.handle.isLargeFile():
                self.residues = sorted(self.residues, 
                                       key=lambda d: int(numonly(d.index)))
        self.indices = sorted(self.indices, key=lambda d: int(numonly(d)))

    def BoundingBoxes(self):
        '''
        Compute the bounding boxes for every residue in this chain.
        '''
        for res in self.GetResidues():
            res._ComputeBoundingBox()

    def pop(self, s):

        """ Pop a residue. """

        res = self.GetResidueByIndex(s)
        if self.structure.handle:
            if not self.structure.handle.isLargeFile():
                self.residues.remove(res)
        self.indices.remove(str(s))
        return res

    def update(self, o):

        """ Add residues from another chain. """

        for res in o: self.AddResidue(res)

    def populate(self):

        """ If a large file, load all residues from memory into this
		object. Otherwise, do nothing. """

        if self.structure:
            handle = self.structure.handle
            if handle and handle.isLargeFile():
                self.residues = self.GetResidues()

    def AsModel(self, name):
        '''
        Return the chain as model

        :param str name: Name of the model
        '''
        model = PDBmodel(name)
        model.structure = self.structure
        model.parent = self.parent
        for r in self.IterResidues():
            r.model = model
            model.AddResidue(r)
        return model

    def FixResidueNumbering(self):
        '''
        If numbering is not consistent in the chain, chainge it to a consecutive one
        '''
        newcha = PDBchain(self.name)
        resind = self.GetIndices()
        for rin in xrange(len(resind)):
            res = self.GetResidueByIndex(resind[rin])
            res.index = str(rin + 1)
            for atom in res.GetAtoms():
                atom.parent = res
            newcha.AddResidue(res)

        return newcha

    def __str__(self):

        get = lambda d: self.GetResidueByIndex(d)
        out = '\n'.join([str(get(x)) for x in self.indices])
        out += '\n%s\n' % (str(PDBterminator(self)))
        return out

    def __getitem__(self, i):

        """ Force the object to load the residue if not present
		in this structure. """

        if self.parent == None:
            parent = -1
        elif isinstance(self.parent, PDBstructure):
            parent = -1
        else:
            parent = self.parent.name
        if not self.structure: self.structure = PDBstructure()
        handle = self.structure.handle
        if str(i) in self.indices:
            if handle and ((handle.isLargeFile()) or (len(self.indices) > \
                                                      len(self.residues))):
                if handle.hasResidue(self.name, str(i), parent):
                    resid = handle.readResidue(self.name, str(i), parent)
                    resid.chain = self.name
                    return resid
                elif self.indices.index(str(i)) < len(self.residues):
                    return self.residues[self.indices.index(str(i))]
                else:
                    raise KeyError(
                        'Could not find residue %s in %s:%s in file or '\
                        'structure.' % (str(i), str(self.name), str(parent)))
            return self.residues[self.indices.index(str(i))]
        else:
            return None

    def __len__(self):
        return len(self.indices)

    def __iter__(self):
        if self.structure == None:
            for ind in self.indices: 
                yield int(ind)
        else:
            if self.structure.islargefile:
                for res in self.indices:
                    yield self.__getitem__(res)            
            else:
                for ind in self.indices: 
                    yield int(ind)                            


class PDBmodel(PDBchain):
    """ A PDB model (a special kind of chain). """

    __slots__ = ('name', 'structure', 'residues', 'indices', 'chains', 
                 'chainNames')

    def __init__(self, name):

        if type(name) != int or name < 0:
            raise ValueError('Model name should be zero or a positive integer.'
                             )
        super(PDBmodel, self).__init__(name)
        self.chains = list()
        self.chainNames = list()
        self.residues = None

    def GetResidues(self):

        """ Acquire all of the residues in this model. Note that if this
		is a large file, this will force the object to load all of these
		residues into memory from the file. """

        this = [self.GetResidueByIndex(i) for i in self.GetIndices()]
        for it in self.chains: this += it.GetResidues()
        return this

    def GetResidueByPosition(self, pos):
        residues = [item for item in self.IterResidues()]
        # residues = self.GetResidues()
        # for chain in self.chains:
        #	residues.extend([res for res in chain.residues])
        return residues[pos]

    def IterResidues(self):

        """ 
        Iteratively yield all of the residues in this model. Note that if this
		is a large file, this will force the object to load all of these
		residues into memory from the file. 
        """

        for i in self.GetIndices(): yield self.GetResidueByIndex(i)
        for i in self.chains:
            for j in i.IterResidues(): yield j

    def GetChains(self):
        return self.chains

    def GetChain(self, name):
        return self.GetChainByName(name)

    def GetChainByName(self, name):
        for x in xrange(len(self.chainNames)):
            if self.chainNames[x] == name: return self.chains[x]

    def GetChainNames(self):
        return self.chainNames

    def NewChain(self, ch):

        """ Create a new PDBchain, add it to this model, and return it. """

        c = PDBchain(ch)
        self.AddChain(c)
        return c

    def AddChain(self, chain):

        """ Add a PDBchain (chain) to this model. """

        chain.structure = self.structure
        chain.parent = self
        self.chains.append(chain)
        self.chainNames.append(chain.name)

    def AddResidue(self, resid):

        """ Add a residue to this model. """

        resid.model = self.name
        if self.residues == None: self.residues = []
        self.residues.append(resid)
        self.indices.append(str(resid.index))

    def CAContactMap(self):
        '''
        Compute the contact map of this model (and any possible chains) using the boundingbox method.
        '''
        contactmap = super(PDBmodel, self).CAContactMap()
        for chain in self.GetChains():  # Deal with children.
            contactmap.extend(chain.CAContactMap())
        return sorted(contactmap, key=lambda e: (e[0], e[1]))

    def BBContactMap(self, thres=4.5):
        '''
        Compute the contact map of this model (and any possible chains) using the boundingbox method.
        '''
        contactmap = super(PDBmodel, self).BBContactMap(thres)
        for chain in self.GetChains():  # Deal with children.
            contactmap.extend(chain.BBContactMap(thres))
        return sorted(contactmap, key=lambda e: (e[0], e[1]))

    def ContactMap(self, thres=4.5, userindices=None):

        """ Compute the contact map of this model (and any possible chains). """

        if userindices:

            contactmap = super(PDBmodel, self).ContactMap(thres=thres, userindices=userindices)
            if not contactmap:
                for chain in self.GetChains():  # Deal with children.
                    contactmap.extend(chain.ContactMap(thres=thres, userindices=userindices))
            return sorted(contactmap, key=lambda e: (e[0], e[1]))

        else:

            contactmap = super(PDBmodel, self).ContactMap(thres)
            if not contactmap:
                for chain in self.GetChains():  # Deal with children.
                    contactmap.extend(chain.ContactMap(thres))
            return sorted(contactmap, key=lambda e: (e[0], e[1]))

    def BoundingBoxes(self):
        '''
        Compute the bounding boxes for every residue in this model.
        '''
        for res in self.GetResidues():
            res._ComputeBoundingBox()

    def populate(self):
        # this needs to be fixed is not for large files!!!! the inheritance 
        # from chain needs to be kept
        if len(self.chains) == 0:
            self.residues = self.GetResidues()
        elif len(self.chains) == 1:
            self.chains[0].populate()
            self.residues = self.chains[0].residues
        else:
            self.residues = {}
            for chain in self.chains:
                chain.populate()
                self.residues[chain.name] = chain.residues

    def __getitem__(self, i):
        handle = self.structure.handle
        if str(i) in self.indices:
            if handle and ((handle.isLargeFile()) or (len(self.indices) > \
                                                      len(self.residues))):
                if handle.hasResidue('', str(i), self.name):
                    resid = handle.readResidue('', str(i), self.name)
                    resid.model = self.name
                    return resid
                elif self.indices.index(str(i)) < len(self.residues):
                    return self.residues[self.indices.index(str(i))]
                else:
                    raise KeyError('Could not find residue %s in MODEL %s in '\
                                   'file or structure.' % (str(i), 
                                                           str(self.name)))
            return self.residues[self.indices.index(str(i))]
        else:
            return None

    def __len__(self):

        l = len(self.indices)
        for c in self.GetChains():
            l += len(c)
        return l

    def __str__(self):
        model = 'MODEL%9s\n' % (str(self.name))
        if not self.GetChainNames():
            #does not have chains
            model += super(PDBmodel, self).__str__()
        else:
            #model has chains
            for ch in self.chains: model += str(ch)
        return model + 'ENDMDL\n'

    def __iter__(self):
        if self.structure == None:
            for res in self.indices: 
                yield int(res)
        else:
            if self.structure.islargefile:
                for res in self.indices:
                    yield self.__getitem__(res)
            else:               
                for res in self.indices: 
                    yield int(res)

class PDBstructure(object):
    """ A PDB protein structure (a collection of chains/models). """

    __slots__ = ('chains', 'orderofchains', 'models', 'orderofmodels', 
                 'remarks', 'filepath', 'organism', 'taxid', 'mutation', 
                 'EC_code', 'DOI', 'PMID', 'contactmap', 'handle', 'ismodel', 
                 'FastaAssociations', 'fastaresiduehomologs', 
                 'FastaChainEquivalence', 'array', 'rmsd_matrix', 
                 'average_rmsd', 'istrajectory','verbose', 'islargefile')

    def __init__(self, filein='', ismodel=False, verbose=False):

        # Attributes.
        self.filepath = filein
        self.handle = None
        self.organism = None
        self.taxid = None
        self.mutation = False
        self.contactmap = list()
        self.orderofchains = list()
        self.orderofmodels = list()
        self.remarks = list()
        self.chains = dict()
        self.models = dict()
        self.ismodel = ismodel
        self.FastaAssociations = {}
        self.fastaresiduehomologs = {}
        self.FastaChainEquivalence = {}
        self.verbose = verbose
        # Read a file?
        if filein != '':
            self.ReadFile(filein)
            self.istrajectory=self.IsTrajectory()

    # Accesssors/Mutators.

    def GetChainNames(self):
        return self.orderofchains

    def GetModelNames(self):
        return self.orderofmodels

    def GetChain(self, ch):

        """ Get a chain by name. """

        if ch in self.chains: return self.chains[ch]
        return None

    def GetModel(self, mod):

        """ Get a model by name. """

        if mod in self.models: return self.models[mod]
        return None

    def NewChain(self, name):

        """ Construct and add a new chain by name to the PDB. Returns the chain. """

        p = PDBchain(name)
        self.AddChain(name, p)
        return p

    def NewModel(self, name):

        """ Construct and add a new model by name to the PDB. Returns the model. """

        p = PDBmodel(name)
        self.AddModel(name, p)
        return p

    def RemoveChain(self, name):

        """ Remove a chain from the structure (by name). Returns the chain. """

        c = self.GetChain(name)
        if not name in self.chains:
            raise KeyError('Chain %s does not exist in structure!' % (name))
        del self.chains[name]
        self.orderofchains.remove(name)
        return c

    def RenameChain(self, name, newname):

        """ Rename a chain in the structure (by name). Returns the chain. """

        c = self.GetChain(name)
        if not name in self.chains:
            raise KeyError('Chain %s does not exist in structure!' % (name))
        elif newname in self.chains:
            raise KeyError('New name %s already exists in structure!' % (name))
        self.chains[newname] = c
        del self.chains[name]
        self.orderofchains[self.orderofchains.index(name)] = newname
        c.name = newname
        return c

    def RemoveModel(self, name):

        """ Remove a model from the structure (by name). Returns the model. """

        m = self.GetModel(name)
        if not name in self.models:
            raise KeyError('MODEL %s does not exist in structure!' % (str(name)))
        del self.models[name]
        self.orderofmodels.remove(name)
        return m

    def RenameModel(self, name, newname):

        """ Renames a model in the structure (by name). Returns the model. """

        m = self.GetModel(name)
        if not name in self.models:
            raise KeyError('MODEL %s does not exist in structure!' % (str(name)))
        self.models[newname] = m
        del self.models[name]
        self.orderofmodels[self.orderofmodels.index(name)] = newname
        m.name = newname
        return m

    def AddChain(self, chainname, chain):

        """ Add a chain as a list of residues to the PDB. """

        # Do chain addition operations.
        if chainname in self.chains:
            raise KeyError('Chain already exists in structure by that name!')
        if type(chain) != PDBchain and type(chain) == PDBmodel:
            # Cast into a PDBchain.
            cast = PDBchain(chainname)
            for i in chain: cast.AddResidue(chain[i])
            chain = cast
        chain.populate()  # Ensure all residues are loaded (if large file).
        chain.name = chainname
        self.chains[chainname] = chain
        self.orderofchains.append(chainname)
        chain.structure = self

    def AddModel(self, modelname, model):

        """ Add a model as a list of residues to the PDB. """

        # Do chain addition operations.
        if model in self.models:
            raise KeyError('Model already exists in structure by that name!')
        # if type(model) != PDBmodel and type(model) == PDBchain:
        # Cast into a PDBmodel.
        # cast = PDBmodel(modelname)
        # for i in model: cast.AddResidue(model[i])
        # model = cast
        cast = PDBmodel(modelname)
        for r in model.GetResidues(): cast.AddResidue(r)
        model.populate()  # Ensure all residues are loaded (if large file).
        model.name = modelname
        self.models[modelname] = model
        self.orderofmodels.append(modelname)
        model.structure = self
        return model

    def AddResidueToChain(self, chain, res):

        """ Add a residue to a chain. Deprecated; use chain class function. """

        ch = self.GetChain(chain)
        if ch != None: ch.AddResidue(res)

    def AddResidueToModel(self, model, res):

        """ Add a residue to a model. Deprecated; use model class function. """

        mo = self.GetModel(model)
        if mo != None: mo.AddResidue(res)

    def AddRemark(self, remark):

        """ Add a remark (note/comment) to the structure/PDB file. """

        if len(remark) == 0:
            self.remarks.append('')
            return
        for it in xrange(0, len(remark), 68):
            self.remarks.append(remark[it:it + 68])

    def GetRemarks(self):

        """ Return all remarks from the PDB as a list of strings. """

        return self.remarks

    def CheckSequential(self):
        '''
        Go in every chain, and re-number the residues in a sequential manner
        '''
        new = PDBstructure(filein='', ismodel=self.ismodel)
        if self.ismodel:
            things = self.orderofmodels
            for mo in things:
                m = self.models[mo]
                m = m.FixResidueNumbering()
                m.structure = new
                m.SortByNumericalIndices()
                new.AddModel(mo, m)
        else:
            things = self.orderofchains
            for ch in things:
                c = self.chains[ch]
                c = c.FixResidueNumbering()
                c.structure = new
                c.SortByNumericalIndices()
                new.AddChain(ch, c)
        return new

    def CheckComplete(self):

        """ For every chain, check to see if every residue is complete (see
        aa_list dictionary). """

        for ch in self.orderofchains:
            for re in self.chains[ch]:
                res = self.chains[ch][re]
                nam = res.name
                if not nam in aa_lists: continue
                check = set(set(aa_lists[nam])).difference(res.atoms)
                if len(check) > 0:
                    return res
        for ch in self.orderofmodels:
            for res in self.models[ch]:
                nam = res.name
                if not nam in aa_lists: continue
                check = set(set(aa_lists[nam])).difference(res.atoms)
                if len(check) > 0:
                    return res
        return True

    # I/O Functionality.

    def _pymol(self):

        # Start pymol if not started.
        import pymol
        pymol.finish_launching()
        return pymol

    def view(self, istrajectory=False):

        """ View the structure in a Pymol window. Requires an installation of Pymol. """

        pym = self._pymol()
        pym.cmd.read_pdbstr(str(self), 'PDBnet')
        if len(self.models) > 0 and not istrajectory:
            pym.cmd.set('all_states', 'on')
        elif istrajectory:
            pass  # SERGIO

    def ViewStructure(self):
        self.view()

    def WriteFile(self, filename):

        """ Write this PDB structure as a single PDB file. """

        fh = open(filename, 'w')
        fh.write(str(self))
        fh.close()

    def read(self, filename):

        """ Alias for ReadFile(). """

        self.ReadFile(filename)

    def ReadFile(self, filename):

        """ Read a PDB file. Populate this PDBstructure. """
        filesize = sizeOfFile(filename)
        
        if filesize >= PDB_LARGE_FILE_SIZE:
            self.islargefile = True
            # Map to memory and read the file.
            p = mappedPDBfile(filename)
            self.handle = p
            self.remarks, self.organism, self.taxid, self.mutation, self.DOI, \
                self.PMID, self.EC_code = p.read()
            largeFile = p.isLargeFile()
    
            # Start creating high-level structures for models/chains.
            if p.hasModels():
                self.ismodel = True
                map(self.NewModel, p.getModelNames())
                for mo, ch, re in p.iterResidueData():
                    model = self.GetModel(mo)
                    if ch != '':
                        if not ch in model.GetChainNames(): 
                            chain = model.NewChain(ch)
                        #chain = model.GetChain(ch)
                        if largeFile:
                            chain.AddIndexOfResidue(re)
                        else:
                            chain.AddResidueByIndex(re)
                    else:
                        if largeFile:
                            model.AddIndexOfResidue(re)
                        else:
                            model.AddResidueByIndex(re)
                for mod in self.GetModelNames():
                    model = self.GetModel(mod)
                    model.structure = self
                    model.SortByNumericalIndices()
                    # Added by Jack
                    model.populate()
                    for chain in model.GetChains():
                        chain.SortByNumericalIndices()
            else:
                self.ismodel = False
                map(self.NewChain, p.getChainNames())
                for _, ch, re in p.iterResidueData():
                    chain = self.GetChain(ch)
                    if largeFile:
                        chain.AddIndexOfResidue(re)
                    else:
                        chain.AddResidueByIndex(re)
                for cha in self.GetChainNames():
                    chain = self.GetChain(cha)
                    chain.structure = self
                    chain.SortByNumericalIndices()
    
            # Remove this component from memory if no longer necessary.
            if not largeFile:
                p.close()
                self.handle = None
                del p
        else:
            self.islargefile = False
            p = PDBfile(filename)
            self.remarks, self.organism, self.taxid, self.mutation, self.DOI,\
                self.PMID, self.EC_code = p.read()
            if p.models != []:
                self.ismodel=True
                map(self.NewModel, p.models)
                for mo, ch, re in p.iterResidueData():
                    model = self.GetModel(mo)
                    if ch != '':
                        if not ch in model.GetChainNames(): 
                            chain = model.NewChain(ch)
                        else:
                            chain=model.GetChain(ch)
                        chain.AddResidue(re)
                    else:
                        model.AddResidue(re)
            else:
                self.ismodel=False
                map(self.NewChain, p.chains)
                for _, ch, re in p.iterResidueData():
                    chain = self.GetChain(ch)
                    chain.AddResidue(re)
                for cha in self.GetChainNames():
                    chain = self.GetChain(cha)
                    chain.structure = self
                    chain.SortByNumericalIndices()                

    def ChainAsFASTA(self, chain):

        """ Return the chain as a FASTA. Deprecated; use chain class function. """

        ch = self.GetChain(chain)
        if ch != None: return ch.AsFASTA()

    def ModelAsFASTA(self, model):

        """ Return the model as a FASTA. Deprecated; use chain class function. """

        mo = self.GetModel(model)
        if mo != None: return mo.AsFASTA()

    # Scoring Functionality

    def _getApropriateResidue(self, chain, position):
        """
        PRIVATE: Given a chain, and the index in the fasta file
        return the apropriate residue. THIS IS TIME CONSUMING
        AND SHOULD BE REVISED
        """
        if self.GetChain(chain):
            ch = self.GetChain(chain)
        else:
            ch = self.GetModel(chain)
        # res = ch.GetResidues()
        if position in ch.fstdict: return ch.fstdict[position]

    def _FastaPair(self, fasta):
        """ Loop over all models or chains and finds association between fasta entries
        and PDBnet entries"""
        # from copy import deepcopy
        associations = {}
        done = []
        found = False
        if isinstance(fasta, str): fasta = FASTAnet.FASTAstructure(fasta, uniqueOnly=False)
        seqs = fasta.orderedSequences
        # seqs2 = deepcopy(seqs)
        ismodel = self.ismodel
        if ismodel:
            things = self.models
        else:
            things = self.chains
        for t in things:
            stseq = things[t].AsFASTA()
            found = False
            f = None
            for f in seqs:
                name = f.name
                if stseq == f.gapLess():
                    if name not in done:
                        done.append(name)
                        associations[t] = name
                        found = True
                        break
        self.FastaAssociations.update(associations)
        if len(associations) != len(seqs):
            fas = FASTAnet.FASTAstructure(uniqueOnly=False)
            for fa in fasta:
                if fa.name in done:
                    fas.addSequence(fa.name, fa.sequence)
        # if len(associations) == len(things): found = True
        # else: found = False

        return found, fasta

    def _iterHomologs(self, fasta, chains=None):
        '''
        This method will create an iterator over homologous residues. If chains are provided, then only
        will return those chains
        :param fasta: fasta structure with the alignment or string
        :param chain: list of chains or models to be dealt with or None if all is required
        '''
        if isinstance(fasta, str): fasta = FASTAnet.FASTAstructure(filein=fasta, uniqueOnly=False)
        if self.ismodel:
            orderofthings = self.orderofmodels
            things = self.models
        else:
            orderofthings = self.orderofchains
            things = self.chains
        if not chains: chains = orderofthings
        # get homologous sites
        if not self.fastaresiduehomologs: homologs = self._FastaPdbMatch(fasta)
        for ch in chains:
            for res in self.fastaresiduehomologs[ch]:
                yield res, ch, things[ch].self.GetResidueByPosition(res)

    def _FastaPdbMatch(self, fasta):
        '''
        This method pairs the homologous residues in a FASTA file
        with their actual residues in a matching .pdb file.
        :param fasta: A FASTAnet object (labblouin)
        '''

        if isinstance(fasta, str):
            fasta = FASTAnet.FASTAstructure(filein=fasta, uniqueOnly=False)
        if self.ismodel:
            orderofthings = self.orderofmodels
            things = self.models
        else:
            orderofthings = self.orderofchains
            things = self.chains

        fastaseqlist = fasta.getSequences()
        # Match each chain/model with its corresponding FASTA sequence in the fasta object
        donelist = []
        matchdict = {}
        # nameeq={}
        for t in orderofthings:
            currthing = things[t]
            currthingseq = ''.join([aa[res.name] for res in
                                    currthing.GetResidues()])
            matchdict[t] = None
            for i in range(len(fastaseqlist)):
                currfastaseq = fastaseqlist[i].gapLess()
                if currthingseq == currfastaseq:
                    if i not in donelist:
                        matchdict[t] = fastaseqlist[i]  # key: thing number // value: fasta sequence
                        # nameeq[t]=fastaseqlist[i].name
                        self.FastaChainEquivalence[t] = fastaseqlist[i].name
                        donelist.append(i)
                        break
                    else:
                        continue
                else:
                    continue

        resdict = {}  # key: chain/model // value: homologous residues

        homologous = fasta.getStrictlyUngappedPositions()  # Find the homologous indices in this FASTA

        # For each chain/model, add its homologous residues to the resdict
        for thing in matchdict:

            if isinstance(matchdict[thing], FASTAnet.FASTAsequence):

                s = matchdict[thing].sequence

                indexdict = {}  # key: index in the ungapped string // value: index in the actual model/chain
                count = 0
                for j in range(len(s)):
                    if s[j].isalpha():
                        indexdict[j] = count
                        count += 1

                finalreslist = [indexdict[index] for index in homologous]

                resdict[thing] = finalreslist

            else:
                print 'mismatch! model/chain {} has no FASTA match'.format(
                    thing)
                resdict[thing] = None

        # Assign the resdict to the class variable
        self.fastaresiduehomologs = resdict
        # return names equivalences
        # return nameeq
        return self.FastaChainEquivalence

    def _asArray(self, fasta, chains=None, CA=False):
        '''
        Transform the structure into a numpy array of shape = (n_chains,
        n_atoms, 3)
        :param fasta: fasta file or :class FASTAnet.FASTAstructure instance
        :param list chains: list of particular chains to be transformed
        :param bool CA: Whether to use alpha carbon or centroid.
        '''
        n_chains = []  # container for atoms in chains
        if CA:
            CA = 'GetCA'
        else:
            CA = 'Centroid'

        ismodel = self.ismodel
        # deal if is chain or model
        if ismodel:
            orderofthings = self.orderofmodels
            things = self.models
        else:
            orderofthings = self.orderofchains
            things = self.chains
        if not self.fastaresiduehomologs: self._FastaPdbMatch(
            fasta)  # check if homologous residues have been identified
        if not chains: chains = orderofthings  # if chains not specified do it 
                                               # for the whole structure
        for ch in chains:
            residues = things[ch].GetResidues()  # get all residues in thing
            # restrict residues to homologous ones
            residues = [residues[x] for x in self.fastaresiduehomologs[ch]]  
            n_atoms = []
            if isinstance(residues, dict):  # if there is chains within models,
                                            # loop over dictionary
                for cha, res in residues.iteritems():
                    ca = getattr(res, CA)()
                    dims = [ca.x, ca.y, ca.z]
                    n_atoms.append(dims)
            elif isinstance(residues, list):  # if residues is a list either 
                                              # self is chain or models have 
                                              # only one chain
                for res in residues:
                    ca = getattr(res, CA)()
                    dims = [ca.x, ca.y, ca.z]
                    n_atoms.append(dims)
            n_chains.append(n_atoms)
        self.array = np.array(n_chains)

        return self.array

    def _iterResidueAssociations(self, fasta, chains=None, fast=None):

        """ 
        PRIVATE: Use an input FASTA file to determine associations between 
        residues as they exist between 2+ chains. Constructs a mapping between 
        chain and sequence and yields all strictly homologous positions as 
        residues, chain-by-chain. DEPRECATED AVOID USING
        """
        found, fasta = self._FastaPair(fasta)

        ismodel = self.ismodel
        if ismodel:
            orderofthings = self.orderofmodels
            things = self.models
        else:
            orderofthings = self.orderofchains
            things = self.chains

        # set pointers if missing
        ks = things.keys()
        if not things[ks[0]].parent:
            for x in ks: things[x].parent = self

        if found:
            seqs = fasta.sequences
        else:
            warnings.warn('Not all structures are used since the sequence ' \
                          'could not be found in the fasta')
            orderofthing = list()
            seqs = dict()
            newthings = {}
            for i in orderofthings:
                for k in self.FastaAssociations:
                    if things[i].name == k:
                        orderofthing.append(i)
                        newthings[i] = things[i]
                        seqs[k] = fasta.sequences[self.FastaAssociations[k]]
            orderofthings = orderofthing
            things = newthings

        if not chains or (len(chains) != len(orderofthings)):
            chains = orderofthings
        chaininds = [orderofthings.index(x) for x in chains]
        if len(chaininds) < 2:
            raise ValueError('Need more than two chains (gave %d).' % (
                len(chaininds)))

        if len(things) != len(seqs):
            raise IOError('FASTA set length not equal to PDB length.')
        posv = fasta.getStrictlyUngappedPositions(chaininds)

        for n in chaininds:  # Go chain by chain and yield information.
            ch = orderofthings[n]
            try:
                S = seqs[self.FastaAssociations[ch]]
            except:
                S = seqs[ch]
            seq = S.sequence
            matches = self.GetFASTAIndices(things[ch], seq)
            residueList = [match for match in matches]
            if residueList[0] is None:
                nm = S.name
                chn = things[ch].name
                raise IOError('No relation found between sequence %s and '%(nm)
                              + 'chain %s.' % (chn))
            residueList = sorted(residueList, key=lambda d: d.fstindex)
            for pos in posv:  # Go position by position.
                yield pos, ch, residueList.pop(0)

    def tmscore(self, fasta, chains=None, native=None, CA=True, fast=None):

        """ Get the TMscore between two chains. Requires a
        FASTA alignment and a value for the length of the
        native structure (e.g., for a pairwise alignment,
        the length of the structure used as a reference
        before alignment was done). The latter is computed
        by assuming the first of both provided chains is the
        native structure; otherwise, uses a provided chain
        name (native input). """

        ismodel = self.ismodel
        if ismodel:
            orderofthings = self.orderofmodels
            things = self.models
        else:
            orderofthings = self.orderofchains
            things = self.chains

        if not chains: chains = orderofthings
        if len(chains) != 2:
            raise ValueError('Need exactly two chains to score.')
        if fast is None: fasta = FASTAnet.FASTAstructure(fasta,
                                                         uniqueOnly=False)
        posvect = fasta.getStrictlyUngappedPositions()
        items = self._iterResidueAssociations(fasta, chains, fast=fast)

        # Get the lengths of the respective chains.
        if native == None:
            leN = len(things[chains[0]])
        else:
            leN = len(things[native])  # Length of the reference structure.
        leT = len(posvect)  # Number of aligned positions.

        # Calculate d_0 for this alignment.
        cuberoot = lambda d: d ** (1. / 3)
        d_0 = 1.24 * cuberoot(leN - 15) - 1.8

        # Get the summation portion of the TMscore.
        sumportion = 0
        posmask = {}
        for pos, _, residue in items:
            if not pos in posmask: posmask[pos] = []
            if CA:
                atom = residue.GetCA()
            else:
                atom = residue.Centroid()
            posmask[pos].append(atom)
        for pos in posvect:  # Assume chA is Native Structure.
            cavect = posmask[pos]
            ca1, ca2 = cavect[0], cavect[1]
            d_i = sqrt((ca1.x - ca2.x) ** 2 +
                       (ca1.y - ca2.y) ** 2 +
                       (ca1.z - ca2.z) ** 2)
            sumdenom = 1 + (d_i / d_0) ** 2
            suminside = 1. / (sumdenom)
            sumportion += suminside

        # Return the TMscore.
        return (1. / leN) * sumportion

    def gdt(self, fasta, chains=None, distcutoffs=[1, 2, 4, 8], CA=True, fast=None):

        """ Get the GDT score between two chains. Requires a
        FASTA alignment. """

        ismodel = self.ismodel

        if ismodel:
            orderofthings = self.orderofmodels
            things = self.models
        else:
            orderofthings = self.orderofchains
            things = self.chains

        if len(distcutoffs) != 4:
            raise ValueError('Need exactly 4 distance cutoff values.')
        if not chains: chains = orderofthings
        if len(chains) != 2:
            raise ValueError('Need exactly two chains to score.')
        items = self._iterResidueAssociations(fasta, chains, fast=fast)

        # Get all relevant atoms.
        posmask = {}
        for pos, _, residue in items:
            if not pos in posmask: posmask[pos] = []
            if CA:
                atom = residue.GetCA()
            else:
                atom = residue.Centroid()
            posmask[pos].append(atom)

        # Get the Euclidean distances between
        # all pairs of aligned positions.
        distances = {}
        for pos in posmask:
            cavect = posmask[pos]
            ca1, ca2 = cavect[0], cavect[1]
            dist = sqrt((ca1.x - ca2.x) ** 2 +
                        (ca1.y - ca2.y) ** 2 +
                        (ca1.z - ca2.z) ** 2)
            distances[pos] = dist

        # Calculate the counts for different cutoffs.
        counts = [0, 0, 0, 0]
        poslen = float(len(distances))
        A, B, C, D = distcutoffs
        for pos in distances:
            dist = distances[pos]
            if dist < A: counts[0] += 1
            if dist < B: counts[1] += 1
            if dist < C: counts[2] += 1
            if dist < D: counts[3] += 1
        GDT_P1 = counts[0] / poslen
        GDT_P2 = counts[1] / poslen
        GDT_P3 = counts[2] / poslen
        GDT_P4 = counts[3] / poslen
        return (GDT_P1 + GDT_P2 + GDT_P3 + GDT_P4) / 4.0

    def rmsd_fast(self, fasta, chains=None, CA='Centroid'):
        '''
        Compute the RMSD using the quaterion-based characteristic
        polynomial approach and multiprocessing.

        :param fasta: fasta file or a :class FASTAnet.FASTAstructure instance
        :type fasta: string or :class FASTAnet.FASTAstructure
        :param chains: chains/models that you wnat to compare.
        :type chains: list. If None (default), all chains/models
        :param CA: string denoting which method to use (case sensitive)
        :type CA: str. 'Centroid' or 'GetCA'
        '''
        if CA == 'GetCA':
            CA = True
        else:
            CA = False
        if chains == None:
            if not self.fastaresiduehomologs: self._FastaPdbMatch(fasta)
            chains = self.fastaresiduehomologs.keys()

        RMSD = {}
        A = self._asArray(fasta, chains=chains, CA=CA)
        combos = list(itertools.combinations(chains, 2))
        rmsds = Parallel(n_jobs=-1)(delayed(_rmsd_qcp)
                                    (A[chains.index(c[0]), :, :],
                                     A[chains.index(c[1]), :, :])
                                    for c in combos)

        for k, v in zip(combos, rmsds): RMSD[k] = v

        self.rmsd_matrix = RMSD

        self.average_rmsd = np.mean(RMSD.values())

        return self.rmsd_matrix

    def rmsd(self, fasta, chains=None, CA=True, fast=None):

        """Get the RMSD between chains. Requires a FASTA alignment."""

        items = self._iterResidueAssociations(fasta, chains, fast=fast)
        # Get all relevant atoms.
        posmask = {}
        for pos, _, residue in items:
            if not pos in posmask: posmask[pos] = []
            residue = self._getApropriateResidue(_, pos)
            if CA:
                atom = residue.GetCA()
            else:
                atom = residue.Centroid()
            posmask[pos].append(atom)
        # Get the Euclidean distance squared between
        # all pairs of aligned positions.
        distsqs = {}
        for pos in posmask:
            cavect = posmask[pos]
            for a in xrange(len(cavect)):
                for b in xrange(a + 1, len(cavect)):
                    if (a, b) not in distsqs: distsqs[(a, b)] = 0
                    ca1, ca2 = cavect[a], cavect[b]
                    distsq = ((ca1.x - ca2.x) ** 2
                              + (ca1.y - ca2.y) ** 2
                              + (ca1.z - ca2.z) ** 2)
                    distsqs[(a, b)] += distsq
        # Calculate the average RMSD.
        rmsd = 0
        for d in distsqs:
            d = distsqs[d]
            r = sqrt((float(d) / len(posmask)))
            rmsd += r
        rmsd /= len(distsqs)

        self.average_rmsd = rmsd
        self.rmsd_matrix = distsqs

        return rmsd

    def rrmsd(self, fasta, chains=None, CA=True, fast=None):

        """ Get the RRMSD between chains. Requires a FASTA alignment.
        See Betancourt & Skolnick, "Universal Similarity Measure for
        Comparison Protein Structures". """
        if CA:
            fastrmsdCA = 'GetCA'
        else:
            fastrmsdCA = 'Centroid'
        if self.ismodel:
            orderofthings = self.orderofmodels
            things = self.models
        else:
            orderofthings = self.orderofchains
            things = self.chains

        if not chains: chains = orderofthings
        if len(chains) != 2:
            raise ValueError('Need exactly two chains to score.')

        R_A = self.RadiusOfGyration([chains[0]])
        R_B = self.RadiusOfGyration([chains[1]])
        avglen = (len(things[chains[0]]) +
                  len(things[chains[1]])) / 2
        L_N, e = avglen - 1, math.e
        c = 0.42 - 0.05 * L_N * e ** (-L_N / 4.7) + 0.63 * e ** (-L_N / 37)
        denom = R_A ** 2 + R_B ** 2 - 2 * c * R_A * R_B
        # alignRMSD = self.rmsd(fasta,chains,CA,fast=fast)
        alignRMSD = self.rmsd(fasta, chains=chains, CA=fastrmsdCA)
        return float(alignRMSD) / sqrt(denom)

    # Other Functionality

    def GetAverage(self, fasta, chains=None, newname=None, closest=None, fast=None):
        """
        Acquire a new chain or model corresponding to an average of all present
        chains or models specified.

        :param closest: Measure to return the closest structure to the average.
        :type clostest: str : rmsd, rrmsd, gdt or tmscore
        """
        current = None
        dummyst = PDBstructure()
        if self.ismodel:
            ismodel = True
            orderofthings = self.orderofmodels
            things = self.models
            if newname == None: newname = 0
            dummyst.AddModel(newname, PDBmodel(newname))
        else:
            ismodel = False
            orderofthings = self.orderofchains
            things = self.chains
            if newname == None: newname = '*'
            dummyst.AddModel(newname, PDBchain(newname))
        if chains == None: chains = orderofthings

        # Make sure all sequences have equal FASTAs.
        if fast is None: fasta = FASTAnet.FASTAstructure(filein=fasta,
                                                         uniqueOnly=False)

        def avgResidue(res_no):

            # Create the dummy residue structure.
            firstch = things[chains[0]]
            firstres = firstch.GetResidues()[res_no]
            fake = PDBresidue(firstres.index, firstres.name)

            # Populate it with atoms.
            atomlist = []
            for a in firstres.GetAtoms():
                atomlist.append(PDBatom(a.serial, a.name, 0, 0, 0,
                                        0, 0, a.symbol, ''))

            # Average all properties.
            count = len(chains)
            for ch in chains:
                res = things[ch].GetResidueByIndex(things[ch].GetIndices()
                                                   [res_no])
                atoms = res.GetAtoms()
                for i in xrange(len(atoms)):
                    atomlist[i].x += atoms[i].x
                    atomlist[i].y += atoms[i].y
                    atomlist[i].z += atoms[i].z
                    atomlist[i].occupancy += atoms[i].occupancy
                    atomlist[i].tempFactor += atoms[i].tempFactor
            for a in atomlist:
                a.x /= float(count)
                a.y /= float(count)
                a.z /= float(count)
                a.occupancy /= float(count)
                a.tempFactor /= float(count)

            for a in atomlist: fake.AddAtom(a)
            return fake

        res_nums = range(len(things[chains[0]]))
        for i in res_nums: dummyst.GetModel(newname).AddResidue(avgResidue(i))
        # dummyst.AddModel(0, output)
        clo = None
        if closest:
            def closer(newname, fasta, things, dummyst, dist=closest,
                       model=ismodel):
                '''
                Get the closest model to the average. dumyst is a PDBstructure
                instance with the mean st.
                '''
                d = 100
                current = None
                mods = [v.name for v in things.itervalues()]
                keys = things.keys()
                # if model: mods = ['Model%d'%(x) for x in mods]
                ugp = fasta.getStrictlyUngappedPositions()
                pair = PDBstructure(ismodel=model)
                if model:
                    mod = dummyst.GetModel(newname)
                    pair.AddModel(newname, mod)
                else:
                    cha = dummyst.GetChain(newname)
                    pair.AddChain(newname, cha)
                outseq = pair.GetModel(newname).AsFASTA()
                if len(outseq) != len(things[things.keys()[0]]):
                    nos = ''
                    for j in xrange(len(als)):
                        for k in ugp:
                            if k == j:
                                nos += outseq[k]
                            else:
                                nos += '-'
                else:
                    nos = outseq
                newfas = FASTAnet.FASTAstructure(uniqueOnly=False)
                newfas.addSequence(newname, nos)
                for i in mods:
                    ke = keys[mods.index(i)]
                    if model:
                        als = fasta.full['Model%d' % (i)]
                    else:
                        als = fasta.full[i]
                    newfas.addSequence(i, als)
                    if model:
                        pair.AddModel(i, things[ke])
                    else:
                        pair.AddChain(i, things[ke])
                    if i != newname:
                        nd = getattr(pair, dist)(newfas, chains=[newname, i],
                                                 fast=True)
                        if dist == 'tmscore': nd = 1 - nd
                        if nd < d:
                            d = nd
                            current = i
                    if model:
                        pair.RemoveModel(i)
                    else:
                        pair.RemoveChain(i)
                    newfas.removeSequence(i)

                if current in things:
                    current = things[current]
                else:
                    current = None
                return dummyst, current

            av, closest = closer(newname, fasta, things, dummyst)
        return av, closest

    def RadiusOfGyration(self, chains=None):

        """
		Acquire the radius of the gyration of the entire, or a portion of, the
		PDB protein molecule.
		"""

        if self.ismodel:
            orderofthings = self.orderofmodels
            things = self.models
        else:
            orderofthings = self.orderofchains
            things = self.chains
        if not chains: chains = orderofthings

        # Get centre of mass of entire molecule.
        molCentroid = [0., 0., 0.]
        numAtoms = 0
        for c in chains:
            for re in things[c]:
                r = things[c][re]
                for at in r.atoms:
                    a = r.atoms[at]
                    molCentroid[0] += a.x
                    molCentroid[1] += a.y
                    molCentroid[2] += a.z
                    numAtoms += 1
        m_x, m_y, m_z = molCentroid
        molCentroid = [m_x / numAtoms, m_y / numAtoms, m_z / numAtoms]

        # Sum all of the distances squared vs. centre of mass.
        distsqr = []
        m_x, m_y, m_z = molCentroid
        for c in chains:
            for re in things[c]:
                r = things[c][re]
                for at in r.atoms:
                    a = r.atoms[at]
                    diff = (a.x - m_x + a.y - m_y + a.z - m_z)
                    distsqr.append(diff ** 2)

        # Return the sum of all of these, divided by number of atoms,
        # square rooted.
        sumdistratio = sum(distsqr) / float(numAtoms)
        return sqrt(sumdistratio)

    def GetAllCentroid(self, chain):

        """ Populates the centroids of all residues. """

        out = []
        if chain:
            if not chain in self.chains:
                raise ValueError('Chain %s does not exist!' % (chain))

            ch = self.chains[chain]

            for res in ch:
                res = ch[res]
                res.Centroid()
                out.append(res.centroid)
        else:
            for res in self.IterAllResidues():
                res.Centroid()
                out.append(res.centroid)

        return out

    def IndexSeq(self, chain, fst):

        """ Store in residues the correct index to the fasta.
		Requires a 1-to-1 correspondance at least a portion
		of the way through. Deprecated; use GetFASTAIndices(). """

        ismodel = self.ismodel

        if ismodel:
            thing = self.GetModel(chain)
        else:
            thing = self.GetChain(chain)
        return self.GetFASTAIndices(thing, fst)

    def GetFASTAIndices(self, thing, fst):

        """
		Given a PDBchain, find 1-to-1 correspondances between it and a FASTA
		sequence.
		"""

        chainseq = thing.AsFASTA()
        ungapped = fst.replace('-', '')
        success = True
        if len(ungapped) != len(chainseq):
            success = False
            yield None
        for i in xrange(0, len(chainseq)):
            # See if there is a failed correspondance.
            if chainseq[i] != ungapped[i]:
                yield None
                success = False
                break
        if success:
            index = -1
            residues = thing.GetResidues()
            for i in residues:
                # i = thing[i]
                index = fst.find(aa[i.name], index + 1)
                i.fstindex = index
                thing.fstdict[index] = i
                yield i

    def IterResiduesFor(self, chains=None, fasta=None):

        """
		Produce an iterator to allow one to iterate over all residues for a
		subset of the structure.
		"""

        if self.ismodel:
            orderofthings = self.orderofmodels
            things = self.models
        else:
            orderofthings = self.orderofchains
            things = self.chains
        if chains == None: chains = orderofthings

        for ch in chains:
            chain = things[ch]
            if fasta:
                for res in chain._iterResidueAssociations(fasta, chains=chain):
                    yield res
            else:
                for res in chain.IterResidues():
                    yield res

    def IterAllResidues(self):

        """
		Produce an iterator to allow one to iterate over all possible residues.
		"""

        for ch in self.chains.keys():
            chain = self.GetChain(ch)
            for res in chain: yield chain[res]
        for mo in self.models.keys():
            model = self.GetModel(mo)
            for res in model: yield res

    def gm(self, fasta, chains=None, atomtype='centroid', typeof='str'):

        """ Acquire Geometric Morphometric data corresponding to the
		(x,y,z) coordinates between all homologous residue positions.
		Requires a FASTA alignment. Options include using alpha-carbon
		positions. By default, uses centroids of residues. Returns a list
		of labels and a list of coordinates as raw GM data. The typeof option
		provides an option for coordinate output; they are returned
		as a semicolon-delimited string (str) or as a numpy 2d array (matrix).
		Atom type can be centroid, CA or all.
		"""

        landmarkinfo = {}

        if self.ismodel:
            orderofthings = self.orderofmodels
            things = self.models
        else:
            orderofthings = self.orderofchains
            things = self.chains
        if chains == None: chains = orderofthings

        labels, coords = [], []  # Output Variables

        # Set up name grabbing from FASTA file for GM labelling.
        # if isinstance(fasta,FASTAnet.FASTAstructure): f = fasta
        # else: f = FASTAnet.FASTAstructure(fasta,uniqueOnly=False)
        # name = lambda d: f.orderedSequences[d].name

        # Acquire all homologous positions as defined in FASTA.
        # items = self._iterResidueAssociations(fasta,chains,f)
        if not self.istrajectory:
            if self.verbose:
                print 'isn\'t a trajectory. doing FastaPDBMatch...'
            items = self._FastaPdbMatch(fasta)
        else:
            if self.verbose:
                line = 'Is a trajectory. Generating data outside of'
                line += 'FastaPDBMatch. generating data...\n'
                line += 'making item dictionary...'
                print line
            items = {model:[int(r.index) for r in self.models[model].residues] 
                     for model in self.models}
            if self.verbose:
                print 'making fastaresiduehomologs...'
            self.fastaresiduehomologs = {chain:[res.index for res in \
                                                self.GetModel(chain).residues]
                                         for chain in orderofthings}
            if self.verbose:
                print 'making FastaChainEquivalence...'
            self.FastaChainEquivalence = {chain: chain for chain in
                                          orderofthings}

        chainsSaw = []
        if typeof == 'matrix':
            row = []
        else:
            row = ''
        # for pos,chain,res in items:
        for chain in orderofthings:
            item = self.FastaChainEquivalence[chain]
            thing = things[chain]
            for r in self.fastaresiduehomologs[chain]:
                if self.istrajectory:
                    res = thing.GetResidueByIndex(r)
                else:
                    res = thing.GetResidueByPosition(r)
                if chain not in chainsSaw:
                    ind = orderofthings.index(chain)
                    nm = item  # name(ind).strip().split(':')[0]
                    if not nm in landmarkinfo:
                        landmarkinfo[nm] = []
                        count = 0
                    # labels.append('>%s:%s' % (nm,nm[:4]))
                    # ^^^ sergio said that this may not be used anymore
                    # (had to do with MATT file formatting)
                    labels.append('>%s' % (nm))
                    if len(chainsSaw) != 0:
                        coords.append(row)
                    chainsSaw.append(chain)
                    if typeof == 'matrix':
                        row = []
                    else:
                        row = ''
                if atomtype == 'CA':
                    atom = res.GetCA()
                    if isinstance(row, str):
                        row += '%f;%f;%f;' % (atom.x,
                                              atom.y,
                                              atom.z)
                    elif isinstance(row, list):
                        row.extend([atom.x, atom.y,
                                    atom.z])
                elif atomtype == 'centroid':
                    atom = res.Centroid()
                    if isinstance(row, str):
                        row += '%f;%f;%f;' % (atom.x,
                                              atom.y,
                                              atom.z)
                    elif isinstance(row, list):
                        row.extend([atom.x, atom.y,
                                    atom.z])
                else:
                    for atom in res.GetAtoms():
                        if isinstance(row, str):
                            row += '%f;%f;%f;' % (atom.x,
                                                  atom.y,
                                                  atom.z)
                        elif isinstance(row, list):
                            row.extend([atom.x, atom.y,
                                        atom.z])

                landmarkinfo[nm].append('%d\t%s\t%s\n' % (count, str(res.index),
                                                          str(res.name)))
                count += 1
        coords.append(row)  # Last chain.

        if typeof == 'matrix': coords = np.array(coords)
        return labels, coords, landmarkinfo

    def makefasta(self):
        '''
        Create a dummy fasta for models
        '''
        names = self.GetModelNames()
        w = open('dummy.fasta', 'w')
        for i in xrange(len(names)):
            w.write('>Model%s\n%s\n' % (names[i], self.ModelAsFASTA(names[i])))
        w.close()

    def WriteGM(self, fasta, gm, chains=None, atomtype='centroid'):

        """
		Write the information present in this PDB between multiple chains as a
		Geometric Morphometric text file. This file will be formatted such that
		individual lines correspond to chains and semi-colons separate the
		(x,y,z) coordinates between all homologous residue positions. Requires
		a FASTA alignment. Options include using alpha-carbon positions.
		By default, uses centroids of residues. atomtype can be centroid, CA
		or all.
		:param fasta: File name of the sequence alignment
		:type fasta: str
		:param gm: File name of the output GM file
		:type gm: str
		:param chains: The chains to be written to gm file
		:type chains: None or str
		:param atomtype: The type of atom information to be extracted.
		:type atomtype: str
		"""
        # if self.ismodel and (fasta == None or fasta == ''):
        #	self.makefasta()
        #	fasta = 'dummy.fasta'
        # Get raw GM data.
        labels, coords, landmarkinfo = self.gm(fasta, chains, atomtype)

        # Start writing.
        with open(gm, 'w') as fgm:
            for ind in xrange(len(labels)):
                fgm.write(labels[ind] + ';' + coords[ind] + '\n')
            # Done.
            # fgm.close()

    def WriteLandmarks2(self, fasta, lm, chains=None):
        """ same as WriteLandmarks, but using _FastaPdbMatch"""
        names = self._FastaPdbMatch(fasta)
        ismodel = self.ismodel
        if ismodel:
            orderofthings = self.orderofmodels
            things = self.models
        else:
            orderofthings = self.orderofchains
            things = self.chains
        if chains == None: chains = orderofthings
        line = ''
        if not chains:
            THINGS = orderofthings
        else:
            THINGS = chains
        for t in THINGS:
            line += '>%s\n' % (names[t])
            th = things[t]
            homologs = self.fastaresiduehomologs[t]
            residues = np.array(th.GetResidues())
            residues = residues[homologs]
            for res in xrange(len(homologs)):  # residues:
                resindex = homologs[res] + 1
                gmindex = res
                resname = residues[res].name
                line += '%d\t%d\t%s\n' % (gmindex, resindex, resname)
        with open(lm, 'w') as F:
            F.write(line)

    def WriteLandmarks(self, fasta, lm, chains=None):

        """ Write the information present in this PDB between multiple
		chains as a landmark text file. This file will be formatted such that
		the file is partitioned in sections starting with chain names and individual
		lines in these correspond to homologous residue positions denoted by
		homologous position, residue number, and residue name tab-delimited. Requires
		a FASTA file. """

        ismodel = self.ismodel

        if ismodel:
            orderofthings = self.orderofmodels
            things = self.models
        else:
            orderofthings = self.orderofchains
            things = self.chains
        if chains == None: chains = orderofthings

        flm = open(lm, 'w')

        # Set up name grabbing from FASTA file for labelling.
        f = FASTAnet.FASTAstructure(fasta, uniqueOnly=False)
        name = lambda d: f.orderedSequences[d].name

        # Acquire all homologous positions as defined in FASTA.
        items = self._iterResidueAssociations(fasta, chains, f)

        # Start writing.
        chainsSaw = []
        ind = -1
        for pos, chain, res in items:
            if chain not in chainsSaw:
                ind = orderofthings.index(chain)
                nm = name(ind).strip().split(':')[0]
                flm.write('>%s\n' % (nm))
                chainsSaw.append(chain)
                ind = 0
            flm.write('%s\t%s\t%s\n' % (ind, res.index, res.name))
            ind += 1

    def FDmatrix(self, fasta, chains=None, scaled=True):
        """
        Compute the form difference matrix (FDM) as explained in Claude 2008.
        It relates to the identification of the most influential residue, with
        respect to the overall shape/structure. If the scaled option is True,
        will return an scaled list (from -1 to 1) of the of lenght equal to the
        number of residues. Otherwise will return the raw FDM, rounded so it
        can be included in a PDB. The scaled version is better for
        vizualization. By default the FDM is computed for all chains, but a
        subset can be passed to the chains option.
        """
        names, data = self.gm(fasta, chains=chains, typeof='matrix')
        # get dimensions
        n_obs = data.shape[0]  # number observations/structures
        n_res = data.shape[1] / 3  # Number of residues
        # Re-sahpe data as a 3D array
        data = data.reshape(n_obs, n_res, 3)
        # compute meanshape
        msh = data.mean(axis=0)
        # compute the form matrix for the mean shape
        FMmsh = cdist(msh, msh)
        FMmsh = np.ma.masked_array(FMmsh, np.isnan(FMmsh))
        # iterate over observations to compute the form difference matrix
        vector = []
        for index in xrange(n_obs):
            ndarray = data[index, :, :]
            FM = cdist(ndarray, ndarray)
            FDM = np.ma.masked_array(FM, np.isnan(FM)) / FMmsh
            utr = FDM[np.triu_indices(FDM.shape[0], k=-1)]
            mfdm = np.absolute(FDM - np.median(utr[~np.isnan(utr)]))
            mdat = np.ma.masked_array(mfdm, np.isnan(mfdm))
            s = np.sum(mdat, axis=1).filled(np.nan)
            vector.append(s)

        # compute FD and scaled if required
        FD = (np.array(vector).sum(axis=0)) / n_obs
        if not scaled:
            return [round(x, 2) for x in (FD)]
        else:
            return [round(x, 3) for x in (2 * (FD - min(FD)) / (max(FD) - min(FD))
                                          - 1)]

    def Map2Protein(self, outname, lis, chain, fasta, betadef=0.00, occdef=0.00):
        """
        Map a list of values (lis), that must have a lenght equal to that of
        the number of residues in the PDB to be mapped (chain). If a list of
        list is provided, the first list will be mapped as the beta factor and
        the second as occupancy
        """
        # Check if more than one thing need to be included
        if len([lis]) == 1:
            lis = [lis]

        # Construct empty PDBstructure.
        dummy = PDBstructure()

        if not self.fastaresiduehomologs:
            _ = self._FastaPdbMatch(fasta)
        # Reset beta and ocuppancy
        if self.models:
            m = self.GetModel(chain)
            newm = dummy.NewModel(m.name)
            for r in m.residues:
                for a in r.GetAtoms():
                    a.tempFactor = betadef
                    if len(lis) > 1:
                        a.occupancy = occdef
                newm.AddResidue(r)
        else:
            m = self.GetChain(chain)
            newm = dummy.NewChain(m.name)
            for r in m.residues:
                for a in r.GetAtoms():
                    a.tempFactor = betadef
                    if len(lis) > 1:
                        a.occupancy = occdef
                newm.AddResidue(r)

            # Acquire all homologous positions as defined in FASTA.
            # items = self._iterResidueAssociations(fasta,None,None)

                # populate new field
                # for pos,chainb,res in items:
            # if not chainb == chain:
            #	continue
        alias = self.fastaresiduehomologs
        for i, pos in enumerate(alias[m.name]):
            res = newm.GetResidueByPosition(pos)
            for a in res.GetAtoms():
                a.tempFactor = lis[0][i]
                if len(lis) > 1:
                    a.occupancy = lis[1][i]

        newm.WriteAsPDB(outname)

    def PopulateBoundingBoxes(self):
        '''
        Populate all of the bounding boxes for every residue in every chain (or
        model) in this protein.
        '''
        if self.ismodel:
            for model in self.models:
                self.models[model].BoundingBoxes()
        else:
            for chain in self.chains:
                self.chains[chain].BoundingBoxes()

    def BBContacts(self, chain=None, thres=4.5):
        '''
        Computes the contacts in this protein based on bounding boxes of the residues.
        :param chain: A list of chain or model names, or a single string (if a chain) or integer (if model).
        None = entire structure.
        :param thres: A threshold for distinguishing contact in Angstroms.
        '''

        from intervaltree import Interval, IntervalTree

        # Populate the bounding boxes for each chain/model
        self.PopulateBoundingBoxes()

        if type(chain) != list and type(chain) != int and type(chain) != str:
            raise ValueError('chain is not an int, string, or list.')
        if type(chain) == list:
            if len(list) == 0:
                raise ValueError('chain is a list of length 0. List must have values in it')

        if self.ismodel:
            if chain == None:  # None = entire structure
                orderofthings = self.orderofmodels  # Use all models
                things = self.models  # All models
            elif type(chain) == list:  # A list of models
                orderofthings = [self.orderofmodels[x] for x in chain]  # List of model names to evaluate
                things = [self.models[name] for name in orderofthings]  # models in the orderofthings list
            elif type(chain) == int:
                things = [self.models[chain]]  # The single model
        else:  # Has chains
            if chain == None:
                orderofthings = self.orderofchains  # Same as above (with chains).......
                things = self.chains
            elif type(chain) == list:
                orderofthings = [self.orderofchains[x] for x in chain]
                things = [self.chains[name] for name in orderofthings]
            elif type(chain) == str:
                things = [self.chains[chain]]

        if len(things) == 1:  # Only one chain/model
            return things[0].BBContactMap(thres=thres)
        else:  # Multiple chains/models
            contactmap = []
            for thing in things:
                contactmap.extend(thing.BBContactMap(thres=thres))
            contactmap = list(set(contactmap))  # remove duplicates...
            return contactmap

    def PCAContacts(self, chain=None):
        '''
        Use PCA to determine contacts for a protein's chain(s).
        :param chain: A list of chain or model names, or a single string (if a chain) or integer (if model).
        None = entire structure.
        :return: contactmap
        '''

        self.PopulateBoundingBoxes()

        if type(chain) != list and type(chain) != int and type(chain) != str:
            raise ValueError('chain is not an int, string, or list.')
        if type(chain) == list:
            if len(list) == 0:
                raise ValueError('chain is a list of length 0. List must have values in it')

        if self.ismodel:
            if chain == None:  # None = entire structure
                orderofthings = self.orderofmodels  # Use all models
                things = self.models  # All models
            elif type(chain) == list:  # A list of models
                orderofthings = [self.orderofmodels[x] for x in chain]  # List of model names to evaluate
                things = [self.models[name] for name in orderofthings]  # models in the orderofthings list
            elif type(chain) == int:
                things = [self.models[chain]]  # The single model
        else:  # Has chains
            if chain == None:
                orderofthings = self.orderofchains  # Same as above (with chains).......
                things = self.chains
            elif type(chain) == list:
                orderofthings = [self.orderofchains[x] for x in chain]
                things = [self.chains[name] for name in orderofthings]
            elif type(chain) == str:
                things = [self.chains[chain]]

        if len(things) == 1:  # Only one chain/model
            return things[0].PCAContactMap()
        else:  # Multiple chains/models
            contactmap = []
            for thing in things:
                contactmap.extend(thing.PCAContactMap())
            contactmap = list(set(contactmap))  # remove duplicates...
            return contactmap

    def CAContacts(self, chain=None):

        if type(chain) != list and type(chain) != int and type(chain) != str:
            raise ValueError('chain is not an int, string, or list.')
        if type(chain) == list:
            if len(list) == 0:
                raise ValueError('chain is a list of length 0. List must have values in it')

        if self.ismodel:
            if chain == None:  # None = entire structure
                orderofthings = self.orderofmodels  # Use all models
                things = self.models  # All models
            elif type(chain) == list:  # A list of models
                orderofthings = [self.orderofmodels[x] for x in chain]  # List of model names to evaluate
                things = [self.models[name] for name in orderofthings]  # models in the orderofthings list
            elif type(chain) == int:
                things = [self.models[chain]]  # The single model
        else:  # Has chains
            if chain == None:
                orderofthings = self.orderofchains  # Same as above (with chains).......
                things = self.chains
            elif type(chain) == list:
                orderofthings = [self.orderofchains[x] for x in chain]
                things = [self.chains[name] for name in orderofthings]
            elif type(chain) == str:
                things = [self.chains[chain]]

        if len(things) == 1:  # Only one chain/model
            return things[0].CAContactMap()
        else:  # Multiple chains/models
            contactmap = []
            for thing in things:
                contactmap.extend(thing.CAContactMap())
            return set(contactmap)

    def Contacts(self, chain=None, thres=4.5, fasta=None):
        '''
        Contacts finds all contacts for a chain/model, a list of chain/models,
        or the entire structure (default). FASTA files can be used to determine
        homologous contacts in an alignment.
        :param chain: the name of a single chain/model, a list of names of
        chains/models, or None, meaning the entire structure.
        :param thres: the threshold distance to determine if two residues
        are in contact.
        :param fasta: A FASTAnet object corresponding to the structure.
        :return: contactmap: a list of tuples of contacts. If a fasta has been
        used, then contactmap will be a dictionary;
        if not, then contactmap is a list.
        if the entire structure is evaluated, contactmap
        will be a set of the list of contacts for all chains/models.
        '''

        if fasta and not self.fastaresiduehomologs:
            self._FastaPdbMatch(fasta)

        if self.ismodel:
            if chain == None:  # None = entire structure
                orderofthings = self.orderofmodels  # Use all models
                things = [self.models[x] for x in self.orderofmodels]
            elif type(chain) == list:  # A list of models
                orderofthings = [self.orderofmodels[x] for x in chain]  # List of model names to evaluate
                things = [self.models[name] for name in orderofthings]  # models in the orderofthings list
            elif type(chain) == int:
                things = [self.models[chain]]  # The single model
        else:  # Has chains
            if chain == None:
                orderofthings = self.orderofchains  # Same as above (with chains).......
                things = [self.chains[x] for x in self.orderofchains]
            elif type(chain) == list:
                orderofthings = [self.orderofchains[x] for x in chain]
                things = [self.chains[name] for name in orderofthings]
            elif type(chain) == str:
                things = [self.chains[chain]]

        if len(things) == 1:  # Only one chain/model
            if self.fastaresiduehomologs:
                return things[0].ContactMap(thres=thres,
                                            userindices=self.fastaresiduehomologs[things[0].name])
            else:
                return sorted(things[0].ContactMap(thres=thres))
        else:  # Multiple chains/models
            if self.fastaresiduehomologs:
                contactmap = {}
                for thing in things:
                    contactmap[thing.name] = thing.ContactMap(thres=thres,
                                                              userindices=self.fastaresiduehomologs[thing.name])
                return contactmap
            else:
                contactmap = []
                for thing in things:
                    contactmap.extend(thing.ContactMap(thres=thres))
                return list(set(sorted(contactmap)))

    def ContactMatrix(self, chain, fname='', thres=4.5):

        '''
		Computes a contact matrix for a single state of a pdb file,
		and saves it as a numpy text file that can be read by other methods.

		:param index: the first actual residue index of the chain to be analyzed (so that array does not go out of bounds)
		:param fname: name of the protein file, so that a properly formatted numpy array text file name can be passed.
		:returns: a numpy matrix of the contacts computed.
		'''

        if self.ismodel:
            thing = self.GetModel(chain)
        else:
            thing = self.GetChain(chain)

        matrix = np.zeros((len(thing), len(thing)))

        cm = thing.ContactMap()

        for pair in cm:
            matrix[pair[0]][pair[1]] = 1
            matrix[pair[1]][pair[0]] = 1

        return matrix

    def FindProteinAttributes(self):

        '''
		Determine if the pdb file consists of chains or models
		to hold different states of the protein.

		:returns: A tuple of the length of the protein chain, the number of states, and the index name of the first residue
		'''

        if len(self.GetModelNames()) != 0:
            plength = self.GetModel(self.GetModelNames()[0]).__len__()  # using models
            pstates = self.GetModelNames()
            # find the index name of the first residue
            # (so that the matrix does not go out of bounds)
            # Because sometimes the first index name is not 0
            start = int(self.GetModel(self.GetModelNames()[0]).GetResidues()[0].index)
            return plength, pstates, start, True
        else:
            plength = self.GetChain(self.GetChainNames()[0]).__len__()  # using chains
            pstates = self.GetChainNames()
            start = int(self.GetChain(self.GetChainNames()[0]).GetResidues()[0].index)
            return plength, pstates, start, False

    def AggregateHomologMatrix(self, matchres, thres=4.5):

        # The matrix will be the size of the matchres list length, which are all the same.
        matrix = np.zeros((len(matchres[matchres.keys()[0]]), len(matchres[matchres.keys()[0]])))

        for i in matchres.keys():

            currhoms = matchres[i]

            mapres = {}

            for j in range(len(currhoms)):
                mapres[currhoms[j]] = j  # key: residue position in actual model
            # value: the index in the matrix

            model = self.models[i]

            cm = model.ContactMap(userindices=currhoms, thres=thres)

            for contact in cm:
                indexA = mapres[contact[0]]  # indexA is the 'x' index in the matrix
                indexB = mapres[contact[1]]  # indexB is the 'y' index in the matrix

                matrix[indexA][indexB] += 1
                matrix[indexB][indexA] += 1

        return matrix

    def AggregateSubsetMatrix(self, subset, thres=4.5):

        # matrix will have the length of the subset list
        matrix = np.zeros((len(subset), len(subset)))

        # matrixpositiondict is used to define which index in the matrix
        # the residue corresponds to.
        matrixpositiondict = {}

        # loop through all residues in the subset
        for i in range(len(subset)):
            matrixpositiondict[subset[i]] = 0

        # loop through all models in the file
        for m in self.orderofmodels:

            model = self.models[m]  # model = the actual model object

            cm = model.ContactMap(userindices=subset, thres=thres)

            for contact in cm:
                indexA = matrixpositiondict[contact[0]]  # indexA is the 'x' index in the matrix
                indexB = matrixpositiondict[contact[1]]  # indexB is the 'y' index in the matrix

                matrix[indexA][indexB] += 1
                matrix[indexB][indexA] += 1

        return matrix

    def AggregateMatrix(self, thres=4.5):

        if self.ismodel:
            pstates = self.orderofmodels
            plength = len(self.models[pstates[0]])
        else:
            pstates = self.orderofchains
            plength = len(self.chains[pstates[0]])
        numstates = len(pstates)

        if self.verbose:
            print "plength: {}".format(plength)
            print "numstates: {}".format(numstates)

        matrix = np.zeros((plength, plength))

        usable_cores = multiprocessing.cpu_count() - 1

        tempstateslist = iter([pstates[x:x + usable_cores] for x in xrange(0, len(pstates), 3)])

        nameslist = []

        parallelizer = Parallel(n_jobs=usable_cores)

        if self.ismodel:  # has models

            for group in tempstateslist:
                g = PDBstructure(ismodel=True)
                for state in group:
                    g.AddModel(state, self.GetModel(state))

                tasks_iterator = (delayed(singlerun)(state, thres, numstates) for state in g)

                matrixlist = parallelizer(tasks_iterator)

                # Sort the matrixlist (it may be out of order)
                matrixlist.sort(key=operator.itemgetter(1))

                for m in matrixlist:
                    matrix = np.add(matrix, m[0])
                    nameslist.append(m[1])

                del g

        else:  # has chains

            for group in tempstateslist:
                g = PDBstructure(ismodel=False)
                for state in group:
                    g.AddChain(state, self.GetChain(state))

                tasks_iterator = (delayed(singlerun)(state, thres, numstates) for state in g)

                matrixlist = parallelizer(tasks_iterator)

                # Sort the matrixlist (it may be out of order)
                matrixlist.sort(key=operator.itemgetter(1))

                for m in matrixlist:
                    matrix = np.add(matrix, m[0])
                    nameslist.append(m[1])

                del g

        translatematrix = np.triu(matrix).T + np.triu(matrix)
        translatematrix[np.diag_indices_from(translatematrix)] /= 2

        return translatematrix, len(pstates)

    def IsWithinResidues(self, contact, length):
        '''
        Check if the contact tuple's values are within the range of residue values.
        :param contact: A contact tuple.
        :return: Boolean if true or false.
        '''

        if contact[0] > length or contact[1] > length:
            return False
        else:
            return True

    def FreqCM(self, thres=4.5, homologset=False, subset=None, fasta=None):

        """
		Computes a frequency matrix.

		:param thres: the threshold distance (in Angstroms) for two residues to be in contact.
		:param homologset: boolean value. set to True if using a homolog structural set,
						   so that only homologous residues across all different structures
						   are evaluated for contact between each other.
		:param subset:	 Initialized to Nonetype, but set equal to a list of residue indices
						 that you want to assess for contact between each other.
						 useful for when you want to only find the contact map
						 between a subset of residues in an MD, for instance.
		"""

        if self.fastaresiduehomologs == None:
            if fasta == None:
                raise Exception('Fasta file not found')
            self._FastaPdbMatch(fasta)

        if homologset:

            results = self.AggregateHomologMatrix(matchres=self.fastaresiduehomologs)
            matrix = results
            numstates = len(self.orderofmodels)

        elif subset:
            if self.verbose:
                print 'in subset freqcm'

            results = self.AggregateSubsetMatrix(subset=subset)
            matrix = results
            numstates = len(self.orderofmodels)

        else:

            # create an aggregate matrix of all of the states' contact matrices
            results = self.AggregateMatrix(thres)
            matrix = results[0]
            numstates = results[1]

        # divide each element in the aggregate contact matrix by the number of states, to get frequency ratio
        if self.verbose:
            print 'div: {}'.format(float(numstates))
        matrix /= (float(numstates))

        return matrix

    def IsTrajectory(self):
        '''
        Check if all models (or chains) in this .pdb file are the same length 
        and sequence.
        '''

        if self.ismodel:
            orderofthings = self.orderofmodels
            ex = self.models[orderofthings[0]]
        else:
            orderofthings = self.orderofchains
            ex = self.chains[orderofthings[0]]

        # loop through the models/chains
        for thingname in orderofthings[1:]:

            if self.ismodel:    other = self.GetModel(thingname)
            else:               other = self.GetChain(thingname)

            # Check for length
            if other.__len__() == ex.__len__():

                # compare the info for each residue
                ex_residues = [res for res in ex.GetResidues()]
                other_residues = [res for res in other.GetResidues()]

                for i in range(len(ex_residues)):

                    ex_res = ex_residues[i]
                    other_res = other_residues[i]

                    if ex_res.name != other_res.name:
                        return False
                    if ex_res.index != other_res.index:
                        return False

            else:
                return False

        return True

    def WriteContacts(self, filename):

        """ Write contact map. """

        fout = open(filename, 'w')
        for a in self.contactmap:
            fout.write('%s\n' % (str(a)))
        fout.close()

    # Internals

    def __iter__(self):

        """ Returns all PDBchains as an iterator. """

        for ch in self.chains: yield self.chains[ch]
        for mo in self.models: yield self.models[mo]

    def __len__(self):

        """ Returns the length of the protein. """

        # For all chains.
        chs = self.chains
        mds = self.models
        return sum([len(chs[x]) for x in chs] + [len(mds[x]) for x in mds])

    def _remarksToString(self):

        remarkstr = ''
        for it in xrange(len(self.remarks)):
            remarkstr += 'REMARK%4s %s\n' % (str(it + 1), self.remarks[it])
        return remarkstr

    def __str__(self):

        """ As a string, outputs structure in the PDB format. """

        out = ''
        out += self._remarksToString()
        for chain in self.orderofchains: out += str(self.GetChain(chain))
        for model in self.orderofmodels: out += str(self.GetModel(model))
        return out


class mappedPDBfile(object):
    __slots__ = ('filePath', 'fileHandle', 'memHandle', 'modelIndices', 
                 'chains', 'residueIndices', 'size')

    def __init__(self, fi):

        self.filePath = fi
        self.fileHandle = open(fi, 'r+b')
        self.size = sizeOfFile(fi)
        self.memHandle = memoryMappedFile(self.fileHandle.fileno(), 0)
        self.residueIndices = dict()
        self.modelIndices = dict()
        self.chains = list()

    def isLargeFile(self):
        """ 
        Return whether or on the PDB file this object represents is incredibly 
        large. 
        """

        return (self.size >= PDB_LARGE_FILE_SIZE)

    def close(self):

        """ Close the file. """

        if self.memHandle:
            self.memHandle.close()
        else:
            self.fileHandle.close()

    def read(self):
        """ 
        Acquire all remarks, and the indices of all models and residues. 
        Returnsmremarks, and biological source information as a tuple (remarks,
        organism, taxid, mutant). 
        """

        # Metadata.
        remarks, organism, taxid, mutant, DOI, PMID, EC = [], '', '', False,\
        '', '', ''

        # Cursors.
        k = 0
        curRes = None
        pos = 0
        line = self.memHandle.readline()
        model = -1

        # Scan the file.
        while line:
            recordType = line[:6].strip()
            if recordType == 'MODEL':
                model = int(line.strip().split()[-1])
                self.modelIndices[model] = pos
            elif recordType == 'ATOM':
                sl = line.split()
                if int(sl[1]) >= 100000:
                    k = 1
                else:
                    k = 0
                chain = line[21 + k].strip()
                if chain not in self.chains:
                    self.chains.append(chain)
                residue_index = line[22 + k:27 + k].strip()
                iCode = line[27].strip()
                residue = residue_index + iCode
                if (curRes == None or residue != curRes):
                    self.residueIndices[(model, chain, residue)] = pos
                    curRes = residue
            elif recordType == 'REMARK':
                remark = line[6:].strip('\n')
                if not remark in remarks: remarks.append(remark)
            elif recordType == 'SOURCE':
                if 'ORGANISM_SCIENTIFIC' in line:
                    organism = '_'.join(line.strip().strip(';').split()[-2:])
                elif 'ORGANISM_TAXID' in line:
                    taxid = line.strip().strip(';').split()[-1]
                elif 'MUTANT' in line or 'MUTATION' in line:
                    mutant = True
            elif 'DOI' in line:
                doi = line.strip().split()[-1]
                if DOI != '':
                    DOI = [DOI]
                    DOI.append(doi)
                else:
                    DOI = doi
            elif 'PMID' in line:
                pmid = line.strip().split()[-1]
                if PMID != '':
                    PMID = [PMID]
                    PMID.append(pmid)
                else:
                    PMID = pmid
            elif 'EC:' in line:
                ec = line.strip().strip(';').split()[-1]
                if EC != '':
                    EC = [EC]
                    EC.append(ec)
                else:
                    EC = ec

            pos = self.memHandle.tell()
            line = self.memHandle.readline()
        self.memHandle.seek(0)

        # Return metadata.
        return (remarks, organism, taxid, mutant, DOI, PMID, EC)

    def hasResidue(self, chain, res, model=-1):
        """ Determine if a residue is present in the PDBfile. """

        return ((int(model), chain, res) in self.residueIndices)

    def readResidue(self, chain, res, model=-1):
        """ 
        Parse a residue from the PDB file and return a PDBresidue. 
        """

        resObj = None
        index = self.residueIndices[(int(model), chain, res)]
        self.memHandle.seek(index)
        line = self.memHandle.readline()
        while line:
            recordType = line[:6].strip()
            if recordType == 'ATOM':
                sl = line.split()
                if int(sl[1]) >= 100000:
                    k = 1
                else:
                    k = 0
                chain = line[21 + k].strip()
                residue = line[22 + k:27 + k].strip() + line[27].strip()
                residue_id = line[17 + k:20 + k].strip()
                if residue == res:
                    if not resObj:
                        resObj = PDBresidue(residue, residue_id)
                    serial = int(line[6 + k:12 + k])
                    atom_name = line[12 + k:16 + k].strip()
                    x = float(line[30 + k:38 + k])
                    y = float(line[38 + k:46 + k])
                    z = float(line[46 + k:54 + k])
                    o = float(line[54 + k:60 + k])
                    b = float(line[60 + k:66 + k])
                    sym = line[76 + k:78 + k].strip()
                    try:
                        charge = line[79 + k:80 + k].strip()
                    except:
                        charge = ' '
                    atom = PDBatom(serial, atom_name, x, y, z, o, b, sym, 
                                   charge)
                    resObj.AddAtom(atom)
                else:
                    break
            line = self.memHandle.readline()
        self.memHandle.seek(0)
        return resObj

    def hasModels(self):
        return (len(self.modelIndices) > 0)

    def getModelNames(self):
        return self.modelIndices.keys()

    def getChainNames(self):
        return self.chains

    def iterResidueData(self):

        """ Yield the model number, chain name, and residue number
        for all residues present in this PDB file, not necessarily
        in order. """

        for mod, ch, res in self.residueIndices.keys():
            yield mod, ch, res

    def getResidueNamesInChain(self, ch):
        res = []
        for c, r in self.residueIndices.keys():
            if c == ch: res.append(r)
        return res

class PDBfile(object):
    """
    PDBfile class without memory mapping
    """
    __slots__ = ('filePath',  'models', 'chains', 'residues', 'data', 'size')

    def __init__(self, fi):

        self.filePath = fi
        self.size = sizeOfFile(fi)
        self.residues = list()
        self.models = list()
        self.chains = list()

    def read(self):
        """ 
        Acquire all remarks, and the indices of all models and residues. 
        Returnsmremarks, and biological source information as a tuple (remarks,
        organism, taxid, mutant). 
        """

        # Metadata.
        remarks, organism, taxid, mutant, DOI, PMID, EC = [], '', '', False,\
            '', '', ''
        
        # Atom data. Will include serial, atomtype, residue, chain, resindex,
        # x, y, z, occupancy, tempfactor, atomname
        self.data=[]
        
        # Cursors.
        k = 0
        curRes = None
        preRes = None
        pos = 0
        model = -1

        # Scan the file.
        with open(self.filePath) as F:
            for line in F:
                recordType = line[:6].strip()
                if recordType == 'MODEL':
                    model = int(line.strip().split()[-1])
                    if not model in self.models:
                        if model != -1:
                            self.models.append(model)
                elif recordType == 'ATOM':
                    sl = line.strip('\n')
                    if len(sl) > 80:
                        k = len(sl) - 80
                    else:
                        k = 0
                    serial= int(line[6:11+k])
                    atom_name = line[12 + k:16 + k].strip()
                    chain = line[21 + k].strip()
                    residue_index = line[22 + k:27 + k].strip()
                    self.residues.append(residue_index)
                    iCode = line[27 + k].strip()
                    residue = residue_index + iCode
                    resName = line[17+k:20+k]
                    curRes = (resName,residue_index)
                    x = float(line[30 + k:38 + k])
                    y = float(line[38 + k:46 + k])
                    z = float(line[46 + k:54 + k])
                    o = float(line[54 + k:60 + k])
                    b = float(line[60 + k:66 + k])
                    sym = line[76 + k:78 + k].strip()
                    try:
                        charge = line[79 + k:80 + k].strip()
                    except:
                        charge = ' '                    
                    atom = PDBatom(serial,atom_name,x,y,z,o,b,sym,charge)
                    if curRes == preRes:
                        Res.AddAtom(atom)
                    elif preRes == None:
                        Res = PDBresidue(index=residue_index, name=resName)
                        Res.AddAtom(atom)
                        preRes = curRes
                    else:
                        if not chain in self.chains:
                            self.chains.append(chain)
                        if not self.data:
                            self.data.append((model, chain, Res))
                        else:
                            if ((model, chain, Res) not in self.data) and Res \
                               not in self.data[-1]:
                                self.data.append((model, chain, Res))
                        Res = PDBresidue(index=residue_index, name=resName)
                        Res.AddAtom(atom)
                        preRes = curRes
                elif recordType == 'TER':
                    if not (model, chain, Res) in self.data and Res not in \
                       self.data[-1]:
                        self.data.append((model, chain, Res))
                elif recordType == 'ENDMDL':
                    if not (model, chain, Res) in self.data and Res not in \
                       self.data[-1]:
                        self.data.append((model, chain, Res))
                elif recordType == 'REMARK':
                    remark = line[6:].strip('\n')
                    if not remark in remarks: remarks.append(remark)
                elif recordType == 'SOURCE':
                    if 'ORGANISM_SCIENTIFIC' in line:
                        organism = '_'.join(
                            line.strip().strip(';').split()[-2:]
                        )
                    elif 'ORGANISM_TAXID' in line:
                        taxid = line.strip().strip(';').split()[-1]
                    elif 'MUTANT' in line or 'MUTATION' in line:
                        mutant = True
                elif 'DOI' in line:
                    doi = line.strip().split()[-1]
                    if DOI != '':
                        DOI = [DOI]
                        DOI.append(doi)
                    else:
                        DOI = doi
                elif 'PMID' in line:
                    pmid = line.strip().split()[-1]
                    if PMID != '':
                        PMID = [PMID]
                        PMID.append(pmid)
                    else:
                        PMID = pmid
                elif 'EC:' in line:
                    ec = line.strip().strip(';').split()[-1]
                    if EC != '':
                        EC = [EC]
                        EC.append(ec)
                    else:
                        EC = ec
        # Return metadata.
        return (remarks, organism, taxid, mutant, DOI, PMID, EC)

    def __iter__(self):
        for e in self.data:
            yield e
    
    def iterResidueData(self):
    
        """ 
        Yield the model number, chain name, residue instance
        for all residues present in this PDB file, not necessarily
        in order. 
        """
        for e in self.data:
            yield e[0], e[1], e[2]
          

# Constants

aa = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
      'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
      'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
      'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
      'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
      'UNK': 'X'}

aa_names = {v: k for k, v in aa.items()}  # Flip above dictionary.

aa_lists = {'ALA': ['N', 'CA', 'C', 'O', 'CB'],
            'CYS': ['N', 'CA', 'C', 'O', 'CB', 'SG'],
            'ASP': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'OD2'],
            'GLU': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'OE2'],
            'PHE': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1',
                    'CE2', 'CZ'],
            'GLY': ['N', 'CA', 'C', 'O'],
            'HIS': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'ND1', 'CD2', 'CE1',
                    'NE2'],
            'ILE': ['N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'CD1'],
            'LYS': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CE', 'NZ'],
            'LEU': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2'],
            'MET': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'SD', 'CE'],
            'ASN': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'ND2'],
            'PRO': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD'],
            'GLN': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'NE2'],
            'ARG': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1',
                    'NH2'],
            'SER': ['N', 'CA', 'C', 'O', 'CB', 'OG'],
            'THR': ['N', 'CA', 'C', 'O', 'CB', 'OG1', 'CG2'],
            'VAL': ['N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2'],
            'TRP': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'NE1',
                    'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
            'TYR': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1',
                    'CE2', 'CZ', 'OH'],
            'UNK': []}

# Debugging

if __name__ == "__main__":
    mystructure = PDBstructure(sys.argv[1])
    print mystructure
    fasta = sys.argv[2]
    if mystructure.ismodel:
        names = mystructure.GetModelNames()
    else:
        names = mystructure.GetChainNames()
    mystructure.WriteGM(fasta, 'test2.gm')
    a, b = mystructure.GetAverage(sys.argv[1].strip('pdb') + 'fasta', 
                                  newname=len(mystructure.models),
                                  closest='tmscore')
    lis = mystructure.FDmatrix(sys.argv[1].strip('pdb') + 'gm')  # for debugging purposes only
    mystructure.Map2Protein('test.pdb', lis, names[0], sys.argv[1])  # for debugging purposes only
    if len(sys.argv) > 2:
        print 'RMSD', mystructure.rmsd(sys.argv[2])
        print 'RRMSD', mystructure.rrmsd(sys.argv[2])
        print 'TMscore', mystructure.tmscore(sys.argv[2])
        print 'GDT', mystructure.gdt(sys.argv[2])
    x = mystructure.GetAllCentroid('A')
    mystructure.Contacts()
    print mystructure.contactmap
