# homology.py
# -------------------------
# May 4, 2012; Alex Safatli
# -------------------------
# Homology Modelling utilities.

import os, math, glob, tempfile, random, subprocess, warnings, sys, cStringIO
import PDBnet as PDBnet
from pfam import extractPDBChain
import __main__
__main__.pymol_argv = ['pymol', '-qc'] # Quiet, no GUI for Pymol.
from IO import *
try:
    import pymol
    from modeller import *
    from modeller.scripts import complete_pdb
except: warnings.warn("Could not locate Pymol or Modeller.",ImportWarning)

# ------------ MANUAL MODELLING -------------- #

class manualModeller:
    '''
    Set up and run Modeller using a given
    folder of PDB files and a target FASTA
    sequence.
    '''
    def __init__(self, outfolder, targetfile, templatefldr):
        self.outFolder = outfolder
        self.targetFile = targetfile
        self.templateFolder = templatefldr
        self.templateFiles = getFilesInFolder(self.templateFolder,'pdb')
        self.__prepare__()
        
    def __buildList__(self):
        '''
        Build a string list of all PDB files.
        '''
        out = ''
        for item in self.templateFiles:
            if (len(out) > 0): out += ", '" + getFileName(item) + "' "
            else: out += "'" + getFileName(item) + "' "
        return out
    
    def __prepare__(self):
        '''
        Prepare input for Modeller.
        '''
        # Make modeller folder.
        if not os.path.isdir(self.outFolder):
            makeFolder(self.outFolder)
            
        # Setup path for alignment PIR file.
        alnfile = os.path.join(self.outFolder,'alignment.seg')
            
        # Setup modeller python script.
        scriptDir = os.path.dirname(os.path.realpath(__file__))
        template = os.path.join(scriptDir,'manual_model')
        modscript = open(template).read()
        modscript = modscript.replace('@@PDBDIR@@', self.templateFolder)
        modscript = modscript.replace('@@ALNFILE@@', 'alignment.seg')
        modscript = modscript.replace('@@KNOWN@@', self.__buildList__())
        modscript = modscript.replace('@@SEQUENCE@@', 'target')
        
        # Write modeller python script.
        modout = os.path.join(self.outFolder,'manual_model.py')
        modfout = open(modout,'w')
        modfout.write(modscript)
        modfout.close()
        
        # Write alignment PIR file.
        p = writeModellerPIR(self.templateFiles,self.targetFile,alnfile)
        if p is None: return
        
    def run(self):
        '''
        Run, call Modeller.
        '''
        orig = os.getcwd()
        try:
            os.chdir(self.outFolder)
            os.system('python manual_model.py')
        finally:
            os.chdir(orig)
            
    def postProcess(self):
        '''
        Make a copy of all found models and rename.
        '''
        # Ideally this should rename on basis of some sort of scoring.
        # However, the actual implementation of this scoring by Biskit
        # still needs to be looked at to keep this consistent. For now,
        # renames sequentially.
        templatePDBs = glob.glob(os.path.join(self.outFolder,'target*.pdb'))
        i = 0
        for pdb in templatePDBs:
            copyFile(pdb,os.path.join(self.outFolder,'model_%.2d.pdb'%(i)))
            i += 1

# ----------------- FASTA -------------------- #

class fasta:
    '''
    Open, read, and organize a FASTA file as an
    object.
    '''
    def __init__(self, fname, removeGaps=True, allUpper=True):
        self.fileName = fname
        self.names = []
        self.sequences = []
        self.read(removeGaps, allUpper)
    def read(self, rgaps, allup):
        # Reads Fasta file.
        fin = open(self.fileName)
        temp = ''
        for line in fin:
            if line[0] == '>':
                self.names.append(line[1:].strip())
                if temp: self.sequences.append(temp)
                temp = ''
            else:
                toAdd = line.strip()
                if rgaps: toAdd = toAdd.replace('-','').replace('.','')
                if allup: toAdd = toAdd.upper()
                temp += toAdd
        self.sequences.append(temp)
        fin.close()

def cleanFastaFile(fname, targetfolder=''):
    '''
    Removes gaps from target FASTA file.
    Allows for BLAST input using this file.
    Creates new files for every sequence entry.
    Returns the filename(s) of the new file(s).
    '''    
    if not os.path.isfile(fname): return None    

    # Collect sequences.

    fastaIn = fasta(fname)

    # Open new file to write to.

    numFiles = 0
    fileNames = []
    seqs = fastaIn.sequences
    names = fastaIn.names
    
    if len(seqs) < 2:
        # Ensures a new file is not made if only
        # one sequence in given FASTA file.
        return [fname]

    for x in xrange(len(seqs)):
        temp = ''.join([s for s in names[x] if \
                        s.isalpha() or s.isdigit()]) # Ensures valid filename.
        temp = temp[0:10] # Ensure less than 10 chars.
        newfilename = checkForFileCollision(\
            os.path.join(targetfolder, temp + '.fasta'),'fasta')
        fileNames.append(newfilename)
        fout = open(newfilename, 'w')
        fout.write('>%s\n'% (names[x]))
        for y in range(0, len(seqs[x]), 80):
            fout.write(seqs[x][y:y+80] + '\n')
        fout.close()

    # Return list of filenames.

    return fileNames

def writeModellerPIR(pdb_files, seq, dest_file):
    out = ""
    for pdb in pdb_files:
        out += ">P1;%s\n" % (getFileName(pdb))
        out += "structure:%s:FIRST    :@:END  :: :  : :  \n*\n" % \
            (getFileName(pdb))
    seqfile = fasta(seq,False,False)
    firstseq = seqfile.sequences[0]
    out += ">P1;target\nsequence: : : : : : : : :\n%s*\n" % (firstseq)
    fout = open(dest_file,'w')
    fout.write(out)
    fout.close()
    if len(out) > 0: return True
    else: return False

def writeAlignmentFromFasta(fasta_file, dest_file):
    '''
    Reads a FASTA file and outputs it to a
    Clustal-like alignment.
    '''    
    fastaIn = fasta(fasta_file, False, False)
    fseqs = fastaIn.sequences
    fnames = fastaIn.names
    fout = open(dest_file, 'w')
    fout.write('align1\n------\n')
    for seq in range(len(fseqs)):
        seqout = fseqs[seq].upper().replace('.','-')
        fout.write('%s\t%s\n' % (fnames[seq], seqout))
    fout.close()
    return dest_file

def writeSOAPSSAlignment(fasta_file, pdb1, pdb2, dest_file):
    '''
    Reads a FASTA file and outputs to a SOAPSS
    alignment file.
    '''
    # Determine what residue numbers the PDBs start at.
    resnum = [getFirstResidueNumber(pdb1), getFirstResidueNumber(pdb2)]
    # Write the alignment file from the FASTA.
    fastaIn = fasta(fasta_file, False, True)
    fseqs = fastaIn.sequences
    fnames = fastaIn.names
    if len(fseqs) != 2:
        raise IndexError('Require FASTA length = 2 for alignment.')
    fout = open(dest_file, 'w')
    seq1 = fseqs[0].replace('.','-')
    seq2 = fseqs[1].replace('.','-')
    seq1out, seq2out = '',''
    if len(seq1) != len(seq2):
        raise IndexError('Need same length sequences in FASTA.')
    for ind in xrange(len(seq1)):
        if seq1[ind] == '-' or seq2[ind] == '-':
            seq1out += seq1[ind].lower()
            seq2out += seq2[ind].lower()
        elif seq1[ind] != '-' and seq2[ind] != '-':
            seq1out += seq1[ind]
            seq2out += seq2[ind]
    fout.write('%d %s\n' % (resnum[0],seq1out))
    fout.write('%d %s\n' % (resnum[1],seq2out))
    fout.close()
    return seq1, seq2

def writeFirstFromFasta(fasta_file, dest_file):
    '''
    Extracts first sequence from a FASTA and
    writes to a new one.
    '''
    fastaIn = fasta(fasta_file, False, False)
    fseq = fastaIn.sequences[0]
    fname = fastaIn.names[0]
    fout = open(dest_file, 'w')
    fout.write('>%s\n'% (fname))
    for y in range(0, len(fseq), 80):
        fout.write(fseq[y:y+80] + '\n')
    fout.close()    

def extractRandomFasta(fasta_file, numRandom, dest_file):
    '''
    Extracts a random number of sequences
    from a FASTA file. Writes to new FASTA.
    '''
    if not os.path.isfile(fasta_file):
        return

    # Collects sequences.

    fastaIn = fasta(fasta_file, False, False)

    # Extract random sequences.

    fseqs = fastaIn.sequences
    fnames = fastaIn.names
    
    numSeqs = len(fseqs)
    seqs = []
    names = []

    indexMap = range(numSeqs)
    random.shuffle(indexMap)
    randMap = indexMap[:numRandom]
    
    for num in randMap:
        seqs.append(fseqs[num])
        names.append(fnames[num])

    # Write to file.

    fout = open(dest_file, 'w')
    for x in xrange(len(seqs)):
        fout.write('>%s\n'% (names[x]))
        for y in range(0, len(seqs[x]), 80):
            fout.write(seqs[x][y:y+80] + '\n')
    fout.close()        

# -------------------------------------------- #

def completePDB(pdbin,pdbout):
    '''
    Given a PDB, clean/complete the PDB using
    Modeller's complete_pdb function.
    '''
    # Squelch stdout.
    save_ = sys.stdout
    sys.stdout = cStringIO.StringIO()
    # Perform PDB completion.
    env = environ()
    env.libs.topology.read('${LIB}/top_heav.lib')
    env.libs.parameters.read('${LIB}/par.lib')    
    m = complete_pdb(env,pdbin)
    m.write(file=pdbout)
    # Reclaim stdout.
    sys.stdout = save_

def getRandomPDBFragment(pdb,k=100):
    '''
    Given a PDB, get a random fragment of length k.
    '''
    pdb = open(pdb)
    dat = [x for x in pdb if x.startswith('ATOM')]
    pdb.close()
    if len(dat) < k: return []
    dat = zip(*[iter(dat)]*k)
    random.shuffle(dat)
    return dat[0]

def getFirstResidueNumber(pdb):
    '''
    Gets the FASTA residue number of the first
    ATOM in a PDB file.
    '''
    pdb = PDBnet.PDBstructure(pdb)
    order = pdb.chains
    firstchain = pdb.chains.keys()[0]
    return int(order[firstchain][0])

def checkForFileCollision(path, ext):
    '''
    Ensures a file is not being overwritten.
    Will recursively keep adding a '0' to end of
    base filename until a unique filename is found.
    Returns the path.
    '''
    if os.path.isfile(path):
        newpath = os.path.join(getFolderName(path),getFileName(path) + '0.' + ext)
        return checkForFileCollision(newpath, ext)
    else: return path
    
def system(instruction):
    '''
    Executes a system command line instruction. Returns
    list of stdout, stderr.
    '''
    sproc = subprocess.Popen(instruction, stdout=subprocess.PIPE, 
                                     stderr=subprocess.PIPE, shell=True)  
    return sproc.communicate()

def cleanModellerFolder(path):
    '''
    Moves intermediate PDBs to another subfolder in the
    Modeller folder output from Biskit.
    '''        
    templatePDBs = glob.glob(os.path.join(path,'target*.pdb'))
    targetFolder = makeFolder(os.path.join(path,'intermediates'))
    for pdb in templatePDBs:
        moveFile(pdb,targetFolder)

def getModellerPDB(workflowfolder, targetfile):
    '''
    Copies the model_00.pdb file from Modeller if present
    in given workflow directory from homologyWorkflow.
    '''
    modeller_folder = os.path.join(workflowfolder, 'modeller')
    if os.path.isdir(modeller_folder):
        modeller_file = os.path.join(modeller_folder, 'model_00.pdb')
        if os.path.isfile(modeller_file):
            copyFile(modeller_file, targetfile)

def compareStructures(foldername1, foldername2, outpath, verbose=False):
    '''
    Given 2 folders with PDB files, calculates the RMSDs
    with all possible pairwise combinations of PDBs; writes
    to file located at path determined by outpath.
    ''' 
    # Running list to ensure no redundant checks.
    crosscheck = []
    # File to write scores to.
    scoreFile = os.path.join(outpath,"structcomparisons.list")
    out = open(scoreFile, "w")

    # Do the actual RMSDs.
    for filesA in glob.glob(os.path.join(foldername1, '*.pdb')):
        for filesB in glob.glob(os.path.join(foldername2, '*.pdb')):
            if (getFileName(filesB), getFileName(filesA)) in crosscheck:
                continue
            if (verbose): print "%s, %s" % (getFileName(filesA), getFileName(filesB))
            line = "%s, %s\t%f\n" % (getFileName(filesA), getFileName(filesB), rmsd(filesA,filesB))
            out.write(line)
            crosscheck.append((getFileName(filesA), getFileName(filesB)))

    # Sorts output RMSDs in descending order.
    sortInstr = 'cat ' + outpath + '/structcomparisons.list | sort -k 3 >> ' \
        + outpath + '/structcomparisons.sorted'
    system(sortInstr)
    return os.path.join(outpath,"structcomparisons.sorted")

def startPymol():
    '''
    Ensures Pymol has been launched.
    '''
    pymol.finish_launching()

def radiusOfGyration(pdbFile,chain,pdb=None):
    '''
    Returns the approximate radius of gyration
    of a given PDB molecule. Deprecated; now in PDBnet.
    '''
    if not pdb: pdb = PDBnet.PDBstructure(pdbFile)
    # Get centre of mass of entire molecule.
    mol_centroid = [0,0,0]
    num_atoms = 0
    for r in pdb.chains[chain]:
        resid = pdb.chains[chain][r]
        for at in resid.atoms:
            a = resid.atoms[at]
            mol_centroid[0] += a.x
            mol_centroid[1] += a.y
            mol_centroid[2] += a.z
            num_atoms += 1
    m_x, m_y, m_z = mol_centroid
    mol_centroid = [m_x/num_atoms,m_y/num_atoms,m_z/num_atoms]
    # Sum all of the distances squared vs. centre of mass.
    distsqr = []
    m_x, m_y, m_z = mol_centroid
    for r in pdb.chains[chain]:
        resid = pdb.chains[chain][r]
        for at in resid.atoms:
            a = resid.atoms[at]
            diff = (a.x-m_x+a.y-m_y+a.z-m_z)
            distsqr.append(diff**2)
    # Return the summation of all distance squared divided by
    # the number of atoms, square rooted.
    sumdistratio = sum(distsqr)/num_atoms
    return math.sqrt(sumdistratio)

def tmscore(alignPDB,alignFASTA):
    
    ''' Returns the TMscore and associated P-value. See Xu &
    Zhang, "How significant is a protein structure similarity with
    TM-score = 0.5?". '''
    
    pdb = PDBnet.PDBstructure(alignPDB)
    if len(pdb.chains) != 2:
        raise IOError('Alignment PDB chain length != 2')
    if not chains: chains = [x for x in pdb.chains][:2]
    tmsc = pdb.tmscore(alignFASTA,chains)
    return tmsc, 1 - math.exp(-math.exp((0.1512-tmsc)/0.0242))

def rrmsd(alignPDB,alignFASTA,rmsd=False,chains=None):
    
    ''' Returns the rRMSD of an aligned pairwise PDB and 
    FASTA file. See Betancourt & Skolnick, "Universal Similarity 
    Measure for Comparison Protein Structures". Deprecated. '''
    
    pdb = PDBnet.PDBstructure(alignPDB)
    if len(pdb.chains) < 2:
        raise IOError('Alignment PDB chain length < 2')
    if not chains:
        chains = [x for x in pdb.chains][:2]

    rr = pdb.rrmsd(alignFASTA,chains)
    if rmsd: return rr, pdb.rmsd(alignFASTA,chains)
    return rr

def rmsd(pdbFile1, pdbFile2):
    '''
    Returns the RMSD of two PDB files.
    Requires: Pymol.
    '''        
    pymol.cmd.load(pdbFile1, 'pdb1')
    pymol.cmd.load(pdbFile2, 'pdb2')
    rmsdOut = pymol.cmd.super('pdb1', 'pdb2')[0]
    pymol.cmd.delete('all')
    return rmsdOut
