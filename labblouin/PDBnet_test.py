import unittest

import labblouin
from labblouin import PDBnet
import math

import copy



class TestPDBatom(unittest.TestCase):

    def setUp(self):

        self.pdb = PDBnet.PDBstructure('PDBnet_test_files/1V6R.pdb')

        self.res = self.pdb.GetModel(self.pdb.orderofmodels[0]).GetResidueByPosition(0)

        self.atom = self.res.atomsOrdered[0]

    def tearDown(self):

        del self.atom
        del self.res
        del self.pdb

    def test_fixname(self):

        dummyname = 'hahaha'

        self.atom.name = dummyname

        self.assertNotEqual(dummyname, self.atom.fixname(),
                         'fixname is not changing an incorrect-lengthed name')

        dummyname = 'H'

        self.atom.name = dummyname

        self.assertEqual(4,len(self.atom.fixname()),
                         'fixname is not padding a single-letter name with 3 spaces')

    def test_Element(self):

        self.atom.name = 'A'

        self.assertEqual('A', self.atom.Element(),
                         'not returning the assigned atom name')

        self.atom.name = 'foofoo'

        self.assertEqual('foofoo', self.atom.Element(),
                     'not returning the assigned atom name')

    def test_ThreeLetters(self):

        self.assertEqual(self.res.name, self.atom.ThreeLetters(),
                        'returning an incorrect parent residue name')

    def test_GetPosition(self):

        pos = (self.atom.x, self.atom.y, self.atom.z)

        self.assertEqual(pos, self.atom.GetPosition(),
                         'not returning the proper x,y,z position from the class attributes')

    def test_DistanceTo(self):

        other = self.res.atomsOrdered[1]

        dis = ((self.atom.x - other.x)**2 + (self.atom.y - other.y)**2 + (self.atom.z - other.z)**2)**0.5

        self.assertEqual(dis, self.atom.DistanceTo(other),
                         'the DistanceTo distance does not equal the Euclidean distance between these atoms')

    def test_Bond(self):

        other = self.res.atomsOrdered[1]

        v = (self.atom.x - other.x, self.atom.y - other.y, self.atom.z - other.z)

        self.assertTupleEqual(v, self.atom.Bond(other),
                         'Bond is not returning the same vector as the sample vector,'
                         'given the same parameters')


class TestPDBresidue(unittest.TestCase):

    def setUp(self):

        self.pdb = PDBnet.PDBstructure('PDBnet_test_files/1V6R.pdb')

        self.res = self.pdb.GetModel(self.pdb.orderofmodels[0]).GetResidueByPosition(0)

    def tearDown(self):

        del self.res
        del self.pdb

    def test_GetAtoms(self):

        atomlist = self.res.atomsOrdered

        self.assertListEqual(atomlist, self.res.GetAtoms(),
                             'the atomsOrdered list and GetAtoms() are not returning the same list')

        for element in self.res.GetAtoms():
            self.assertIs(type(element), labblouin.PDBnet.PDBatom,
                         'the list contains elements of type {}, not type(PDBnet.PDBatom) <-- they should be this'
                         .format(type(element)))

    def test_AddAtom(self):

        x = copy.deepcopy(self.res)
        a = copy.deepcopy(self.res.atomsOrdered[1])

        # make a 666 positioned atom.
        a.x = 6
        a.y = 6
        a.z = 6
        a.name = 'NEW'

        x.atoms[a.name.strip()] = a
        x.atomsOrdered.append(a)
        a.parent = x

        self.res.AddAtom(a)

        '''
        self.addTypeEqualityFunc(PDBnet.PDBresidue, self.assertEqual)

        self.assertEqual(x, self.res)
        '''

        self.assertEqual(x.atoms.keys(), self.res.atoms.keys(),
                         'the atoms dictiorary\'s keys are not the same values')

        for key in self.res.atoms.keys():

            self.assertEqual(x.atoms[key].name, self.res.atoms[key].name,
                             'added atom\'s name is not equal in the atoms dict')
            self.assertEqual(x.atoms[key].serial, self.res.atoms[key].serial,
                             'added atom\'s serial is not equal in the atoms dict')
            self.assertEqual(x.atoms[key].x, self.res.atoms[key].x,
                             'added atom\'s x is not equal in the atoms dict')
            self.assertEqual(x.atoms[key].y, self.res.atoms[key].y,
                             'added atom\'s y is not equal in the atoms dict')
            self.assertEqual(x.atoms[key].z, self.res.atoms[key].z,
                             'added atom\'s z is not equal in the atoms dict')
            self.assertEqual(x.atoms[key].occupancy, self.res.atoms[key].occupancy,
                             'added atom\'s occupancy is not equal in the atoms dict')
            self.assertEqual(x.atoms[key].tempFactor, self.res.atoms[key].tempFactor,
                             'added atom\'s tempFactor is not equal in the atoms dict')
            self.assertEqual(x.atoms[key].charge, self.res.atoms[key].charge,
                             'added atom\'s charge is not equal in the atoms dict')
            self.assertEqual(x.atoms[key].symbol, self.res.atoms[key].symbol,
                             'added atom\'s symbol is not equal in the atoms dict')

        self.assertEqual(len(x.atomsOrdered), len(self.res.atomsOrdered),
                         'the atomsOrdered list is an improper length')

        for i in range(len(self.res.atomsOrdered)):

            self.assertEqual(x.atomsOrdered[i].name, self.res.atomsOrdered[i].name,
                             'added atom\'s name is not equal in the atomsOrdered list')
            self.assertEqual(x.atomsOrdered[i].serial, self.res.atomsOrdered[i].serial,
                             'added atom\'s serial is not equal in the atomsOrdered list')
            self.assertEqual(x.atomsOrdered[i].x, self.res.atomsOrdered[i].x,
                             'added atom\'s x is not equal in the atomsOrdered list')
            self.assertEqual(x.atomsOrdered[i].y, self.res.atomsOrdered[i].y,
                             'added atom\'s y is not equal in the atomsOrdered list')
            self.assertEqual(x.atomsOrdered[i].z, self.res.atomsOrdered[i].z,
                             'added atom\'s z is not equal in the atomsOrdered list')
            self.assertEqual(x.atomsOrdered[i].occupancy, self.res.atomsOrdered[i].occupancy,
                             'added atom\'s occupancy is not equal in the atomsOrdered list')
            self.assertEqual(x.atomsOrdered[i].tempFactor, self.res.atomsOrdered[i].tempFactor,
                             'added atom\'s tempFactor is not equal in the atomsOrdered list')
            self.assertEqual(x.atomsOrdered[i].charge, self.res.atomsOrdered[i].charge,
                             'added atom\'s charge is not equal in the atomsOrdered list')
            self.assertEqual(x.atomsOrdered[i].symbol, self.res.atomsOrdered[i].symbol,
                             'added atom\'s symbol is not equal in the atomsOrdered list')

    def test_GetCA(self):

        the_return = self.res.GetCA()

        self.assertIs(type(the_return), labblouin.PDBnet.PDBatom,
                              'GetCA is not returning a PDBnet.PDBatom object')

        self.assertEqual('CA', self.res.GetCA().name,
                         'GetCA is not returning the alpha carbon (the name isn\'t \'CA\'')

        # testing for when an atom is deleted.
        del self.res.atoms['CA']
        for a in self.res.atomsOrdered:
            if a.name == 'CA':
                self.res.atomsOrdered.remove(a)

        with self.assertRaises(LookupError):
            self.res.GetCA()

    def test_Centroid(self):

        r = copy.deepcopy(self.res)

        x = 0
        y = 0
        z = 0

        for atom in r.atoms:
            if atom in ['C','N','O']:
                continue
            x += r.atoms[atom].x
            y += r.atoms[atom].y
            z += r.atoms[atom].z

        div = float(len(r.atoms) - 3)

        x /= div
        y /= div
        z /= div

        r.centroid = PDBnet.PDBatom(None, 'Centroid', x, y, z, 0, 0, '', '')
        r.centroid.parent = self

        self.assertIs(type(self.res.Centroid()), labblouin.PDBnet.PDBatom,
                      'Centroid is returning an object of type {}, not type(PDBnet.PDBatom)'
                      .format(type(self.res.Centroid)))

        self.assertEqual(r.Centroid().x, self.res.Centroid().x,
                         'Centroid is returning an atom with an incorrect x value')

        self.assertEqual(r.Centroid().y, self.res.Centroid().y,
                         'Centroid is returning an atom with an incorrect y value')

        self.assertEqual(r.Centroid().z, self.res.Centroid().z,
                         'Centroid is returning an atom with an incorrect z value')

        # change the backbone atom positions, and see if it affects the centroid position.
        # backbone atoms SHOULD NOT participate in the centroid position calculation,
        # so this should have no effect on the final outcome.

        for key in self.res.atoms:
            if key in ['C','N','O']:
                self.res.atoms[key].x = 0
                self.res.atoms[key].y = 0
                self.res.atoms[key].z = 0

        self.assertEqual(r.Centroid().x, self.res.Centroid().x,
                         'After changing backbone atom positions, Centroid is returning a different'
                         'centroid atom x position (it should not)')
        self.assertEqual(r.Centroid().y, self.res.Centroid().y,
                         'After changing backbone atom positions, Centroid is returning a different'
                         'centroid atom y position (it should not)')
        self.assertEqual(r.Centroid().z, self.res.Centroid().z,
                         'After changing backbone atom positions, Centroid is returning a different'
                         'centroid atom z position (it should not)')

        # change the side-chain atom positions, and see if it affects the centroid position.
        # side-chain atoms SHOULD participate in the centroid position calculation,
        # so this should have an effect on the final outcome.

        for key in self.res.atoms:
            if key not in ['C', 'N', 'O']:
                self.res.atoms[key].x = 0
                self.res.atoms[key].y = 0
                self.res.atoms[key].z = 0

        self.assertNotEqual(r.Centroid().x, self.res.Centroid().x,
                            'After changing side-chain atom positions, Centroid is still returning'
                            'the same centroid atom x position (it should not)')

        self.assertNotEqual(r.Centroid().y, self.res.Centroid().y,
                            'After changing side-chain atom positions, Centroid is still returning'
                            'the same centroid atom y position (it should not)')

        self.assertNotEqual(r.Centroid().z, self.res.Centroid().z,
                            'After changing side-chain atom positions, Centroid is still returning'
                            'the same centroid atom z position (it should not)')

    def test_BondAngle(self):

        x = (0.0,1.0,2.0)
        y = (3.0,4.0,5.0)

        dp = x[0]*y[0] + x[1]*y[1] + x[2]*y[2]

        l1 = ((x[0]-0.0)**2 + (x[1]-0.0)**2 + (x[2]-0.0)**2)**0.5
        l2 = ((y[0]-0.0)**2 + (y[1]-0.0)**2 + (y[2]-0.0)**2)**0.5

        angle = math.acos(dp / (l1*l2)) * 180 / math.pi

        function_result = self.res.BondAngle((0.0,1.0,2.0),(3.0,4.0,5.0))

        self.assertEqual(angle, function_result,
                         ' BondAngle is returning an improper value. Should be {}, returning {}'
                         .format(angle, function_result))

        self.assertIs(type(function_result), float,
                      'BondAngle returns a {}, not a float'.format(type(function_result)))

    def test_HydrogenBond(self):
        self.fail('waiting for Michelle\'s commit')

    def test_InContactWith(self):

        firstmodel = self.pdb.orderofmodels[0]

        for i in range(len(self.pdb.GetModel(firstmodel).residues)):

            other = self.pdb.GetModel(firstmodel).residues[i]

            # assuming the DistanceTo test passed....
            d = self.res.atoms['CA'].DistanceTo(other.atoms['CA'])

            if d <= 4.5:
                self.assertEqual(True, self.res.InContactWith(other),
                                 'InContactWith is returning a False value when the CA distances are less'
                                 'than 4.5 (should be true)')
            elif d < 12:
                for atomA in self.res.atoms:
                    if atomA in ['C', 'N', 'O']:
                        continue
                    atomA = self.res.atoms[atomA]

                    for atomB in other.atoms:
                        if atomB in ['C', 'N', 'O']:
                            continue
                        atomB = other.atoms[atomB]

                        if atomA.DistanceTo(atomB) <= 4.5:
                            self.assertEqual(True, self.res.InContactWith(other),
                                             'InContactWith is returning a False values when an atom-pair is'
                                             ' less than 4.5 apart (should be true)')
            else:  # They're not in contact
                self.assertEqual(False, self.res.InContactWith(other),
                                 'InContactWith is returning a True value when the CA distances are >= 12'
                                 ' apart (should be False)')


    def test__ComputeBoundingBox(self):
        self.fail()


class TestPDBchain(unittest.TestCase):

    def setUp(self):

        self.pdb = PDBnet.PDBstructure('PDBnet_test_files/1V6R.pdb')

        self.chain = self.pdb.GetModel(self.pdb.orderofmodels[0]).chains[0]

    def tearDown(self):

        del self.pdb
        del self.chain

    def test_GetResidues(self):

        for r in self.chain.GetResidues():

            print type(r)

            self.assertIs(type(r), labblouin.PDBnet.PDBresidue,
                          'GetResidues is returning a list with an element of type {}: '
                          'they all need to be type(PDBnet.PDBresidue)'
                          .format(type(r)))

    def test_IterResidues(self):
        self.fail()

    def test_GetResidueByIndex(self):
        self.fail()

    def test_GetResidueByPosition(self):
        self.fail()

    def test_GetAtoms(self):
        self.fail()

    def test_GetIndices(self):
        self.fail()

    def test_AddIndexOfResidue(self):
        self.fail()

    def test_AddResidue(self):
        self.fail()

    def test_AddResidueByIndex(self):
        self.fail()

    def test_RemoveResidue(self):
        self.fail()

    def test_PCAContactMap(self):
        self.fail()

    def test_BBContactMap(self):
        self.fail()

    def test_CAContactMap(self):
        self.fail()

    def test_ContactMap(self):
        self.fail()

    def test_GetPrimaryPropertiesFromBioPython(self):
        self.fail()

    def test_AsFASTA(self):
        self.fail()

    def test_WriteAsPDB(self):
        self.fail()

    def test_SortByNumericalIndices(self):
        self.fail()

    def test_BoundingBoxes(self):
        self.fail()

    def test_pop(self):
        self.fail()

    def test_update(self):
        self.fail()

    def test_populate(self):
        self.fail()

    def test_AsModel(self):
        self.fail()

    def test_FixResidueNumbering(self):
        self.fail()


class TestPDBmodel(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_GetResidues(self):
        self.fail()

    def test_GetResidueByPosition(self):
        self.fail()

    def test_IterResidues(self):
        self.fail()

    def test_GetChains(self):
        self.fail()

    def test_GetChain(self):
        self.fail()

    def test_GetChainByName(self):
        self.fail()

    def test_GetChainNames(self):
        self.fail()

    def test_NewChain(self):
        self.fail()

    def test_AddChain(self):
        self.fail()

    def test_AddResidue(self):
        self.fail()

    def test_CAContactMap(self):
        self.fail()

    def test_BBContactMap(self):
        self.fail()

    def test_ContactMap(self):
        self.fail()

    def test_BoundingBoxes(self):
        self.fail()

    def test_populate(self):
        self.fail()


class TestPDBstructure(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_GetChainNames(self):
        self.fail()

    def test_GetModelNames(self):
        self.fail()

    def test_GetChain(self):
        self.fail()

    def test_GetModel(self):
        self.fail()

    def test_NewChain(self):
        self.fail()

    def test_NewModel(self):
        self.fail()

    def test_RemoveChain(self):
        self.fail()

    def test_RenameChain(self):
        self.fail()

    def test_RemoveModel(self):
        self.fail()

    def test_RenameModel(self):
        self.fail()

    def test_AddChain(self):
        self.fail()

    def test_AddModel(self):
        self.fail()

    def test_AddResidueToChain(self):
        self.fail()

    def test_AddResidueToModel(self):
        self.fail()

    def test_AddRemark(self):
        self.fail()

    def test_GetRemarks(self):
        self.fail()

    def test_CheckSequential(self):
        self.fail()

    def test_CheckComplete(self):
        self.fail()

    def test__pymol(self):
        self.fail()

    def test_view(self):
        self.fail()

    def test_ViewStructure(self):
        self.fail()

    def test_WriteFile(self):
        self.fail()

    def test_read(self):
        self.fail()

    def test_ReadFile(self):
        self.fail()

    def test_ChainAsFASTA(self):
        self.fail()

    def test_ModelAsFASTA(self):
        self.fail()

    def test__getApropriateResidue(self):
        self.fail()

    def test__FastaPair(self):
        self.fail()

    def test__iterHomologs(self):
        self.fail()

    def test__FastaPdbMatch(self):
        self.fail()

    def test__asArray(self):
        self.fail()

    def test__iterResidueAssociations(self):
        self.fail()

    def test_tmscore(self):
        self.fail()

    def test_gdt(self):
        self.fail()

    def test_rmsd_fast(self):
        self.fail()

    def test_rmsd(self):
        self.fail()

    def test_rrmsd(self):
        self.fail()

    def test_GetAverage(self):
        self.fail()

    def test_RadiusOfGyration(self):
        self.fail()

    def test_GetAllCentroid(self):
        self.fail()

    def test_IndexSeq(self):
        self.fail()

    def test_GetFASTAIndices(self):
        self.fail()

    def test_IterResiduesFor(self):
        self.fail()

    def test_IterAllResidues(self):
        self.fail()

    def test_gm(self):
        self.fail()

    def test_makefasta(self):
        self.fail()

    def test_WriteGM(self):
        self.fail()

    def test_WriteLandmarks2(self):
        self.fail()

    def test_WriteLandmarks(self):
        self.fail()

    def test_FDmatrix(self):
        self.fail()

    def test_Map2Protein(self):
        self.fail()

    def test_PopulateBoundingBoxes(self):
        self.fail()

    def test_BBContacts(self):
        self.fail()

    def test_PCAContacts(self):
        self.fail()

    def test_CAContacts(self):
        self.fail()

    def test_Contacts(self):
        self.fail()

    def test_ContactMatrix(self):
        self.fail()

    def test_FindProteinAttributes(self):
        self.fail()

    def test_AggregateHomologMatrix(self):
        self.fail()

    def test_AggregateSubsetMatrix(self):
        self.fail()

    def test_AggregateMatrix(self):
        self.fail()

    def test_IsWithinResidues(self):
        self.fail()

    def test_FreqCM(self):
        self.fail()

    def test_IsTrajectory(self):
        self.fail()

    def test_WriteContacts(self):
        self.fail()

    def test__remarksToString(self):
        self.fail()


class TestPDBfile(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_isLargeFile(self):
        self.fail()

    def test_close(self):
        self.fail()

    def test_read(self):
        self.fail()

    def test_hasResidue(self):
        self.fail()

    def test_readResidue(self):
        self.fail()

    def test_hasModels(self):
        self.fail()

    def test_getModelNames(self):
        self.fail()

    def test_getChainNames(self):
        self.fail()

    def test_iterResidueData(self):
        self.fail()

    def test_getResidueNamesInChain(self):
        self.fail()


# -----------------------------------------------

if __name__ == '__main__':
    unittest.main()