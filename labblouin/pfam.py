# pfam.py
# -------------------------
# May 14, 2012; Alex Safatli
# -------------------------
# Operations for PFAM.

import labblouin.PDBnet as PDBnet
import xml.dom.minidom
import urllib2, urllib, os, time, gzip, math
try: import texttable
except: pass

# Construct an object representing a PFAM flat file (Stockholm file format).

class pfamFile:
    
    def __init__(self, filepath):
        self.filePath = filepath
        self.clusterResults = None
        self.sequences = {}
        self.clusterSequences = {}
        self.structures = {}
        self.accessionNums = {}
        self.numLines = 0
    
    def parse(self):
        # Does the actual parsing of the file.
        fin = open(self.filePath, 'r')
        line = 'null'
        lineno = 1
        while (line != ''):
            line = self.__parseline(fin.readline())
            if lineno == 1 and line != "STOCKHOLM 1.0":
                print 'ERROR: File given is not a Stockholm file.'
                return
            lineno += 1
        self.numLines = lineno
        
    def __parseline(self, line):
        # Parses an individual line.
        extr = '' # Extracted data.
        
        if line.startswith("#"):
            # Is an ANNOTATION.
            extr = line[1:].strip()
            if extr == "STOCKHOLM 1.0": pass
            elif extr.startswith("=GS"):
                ann = extr[3:].strip().split()
                if ann[1] == 'AC':
                    # Is an accession number for a given sequence.
                    self.accessionNums[ann[0]] = ann[2]
                elif ann[1].startswith('DR'):
                    # Is structure information (extracts: PDB ID, chain).
                    self.structures[ann[0]] = (ann[3],ann[4].split(';')[0])
        else:
            # Is not an annotation -- is an ALIGNMENT or empty string.
            extr = line
            if not (len(line) == 0): # If not an empty string.
                seq = line.split()
                if len(seq) == 2: self.sequences[seq[0]] = seq[1]
                
        return extr

    def cluster(self):
        # Goes through all sequences and 
        # clusters them (i.e. sees which are
        # identical). Useful for full
        # alignments. Assumes ordered by tree.
        seq_list = {}
        curr_seq = ''
        curr_num = len(self.sequences)
        for key in self.sequences:
            if not curr_seq == self.sequences[key]:
                seq_list[key] = self.sequences[key]
                curr_seq = self.sequences[key]                
        self.sequences = seq_list
        new_num = len(seq_list)
        self.clusterResults = (curr_num, new_num)
        
# Other Functions

def extractPDBSequences(pdb_file):
    '''
    Extracts the Amino Acid sequences for a protein structure
    defined by a PDB file and returns as a dictionary.
    '''
    out_seq_list = {}
    
    pdb = PDBnet.PDBstructure(pdb_file)
    for ch in pdb.chains:
        # For every chain in the PDB file
        # get the FASTA sequence.
        seq = pdb.ChainAsFASTA(ch)
        out_seq_list[ch] = seq
        
    return out_seq_list

def extractPDBChain(pdb_file, ch, dest_file):
    '''
    Extracts a PDB file from another that has only
    information for a single chain.
    '''
    fin = open(pdb_file,'r')
    out = []
    for line in fin:
        if not line.startswith('ATOM'): continue
        sl = line.split()
        if int(sl[1]) >= 100000: k = 1
        else: k = 0
        chain = line[21+k].strip()
        if chain == ch: out.append(line)
    fin.close()
    fout = open(dest_file, 'w')
    for x in out: fout.write(x)
    fout.close()    
    return dest_file

def writeSequencesToFile(seq_list, targetfolder, fasta=True):
    '''
    Takes in a dictionary of sequences and writes them
    to seperate .seq or .fasta files.
    '''    
    out_file_list = []

    for ch, seq in seq_list.iteritems():
        if (fasta): ext = 'fasta'
        else: ext = 'sequence'
        path = os.path.join(targetfolder, ch + '.' + ext)
        fout = open(path, "w")
        if (fasta): fout.write('>%s\n' % ch)
        for y in range(0, len(seq), 80): fout.write(seq[y:y+80] + '\n')
        fout.close()
        out_file_list.append(path)
        
    return out_file_list

def writeSequenceToFasta(seq, seq_id, targetpath):
    '''
    Creates FASTA file from given sequence string.
    '''
    fout = open(targetpath, 'w')
    fout.write('>%s\n' % (seq_id))
    for y in range(0, len(seq), 80):
        fout.write(seq[y:y+80] + '\n')
    fout.close()
    
def writeSequencesToFasta(seq_list, dest_file):
    '''
    Creates single FASTA file from given sequence tuple list
    where first item is sequence ID, second item is the sequence.
    '''
    fout = open(dest_file, "w")
    for s in seq_list:
        seq = s[1]
        ch = s[0]
        fout.write('>%s\n' % ch)
        for y in range(0, len(seq), 80):
            fout.write(seq[y:y+80] + '\n')
    fout.close()

def doSequenceSearch(seq_file):
    '''
    Performs the PFAM search for a sequence file.
    '''       
    # Set up curl request.
    url = 'http://pfam.sanger.ac.uk/search/sequence'
    values = {'seq': open(seq_file).read(), 'evalue': '1.0', 'ga':'0', 'output':'xml'}
    data = urllib.urlencode(values)
    req = urllib2.Request(url,data)
    # Retrieve search results.
    response = urllib2.urlopen(req).read()
    dom = xml.dom.minidom.parseString(response)
    result_url = dom.getElementsByTagName('result_url')[0].toxml()
    # Find URL for results.
    result_url = result_url.replace('<result_url>','').replace('</result_url>','')
    # Get results.
    results = ''
    while results == '': # Ensures job is done server-side.
        result_out = urllib2.Request(result_url)
        results = urllib2.urlopen(result_out).read()
    return results

def processSequenceSearch(search_results):
    '''
    Takes in a XML result returned from a sequence
    search and parses it in order to determine what
    the most significant results are for related families.
    Returns a list of the family PFAM IDs.
    '''        
    sigresults = []
    dom = xml.dom.minidom.parseString(search_results)
    matches = dom.getElementsByTagName('location')
    for match in matches:
        if match.getAttribute('significant') == '1':
            sigresult = match.parentNode.getAttribute('id')
            sigresults.append(sigresult)
    return sigresults

def downloadFamilySequences(pfam_family_id, dest_folder, atype='seed'):
    '''
    Takes in PFAM family ID and acquires the gzipped
    flat file filled with seed or full sequence allignments.
    '''     
    dest = os.path.join(dest_folder, '%s-pfam-%s.gz' % (atype,pfam_family_id))
    if os.path.isfile(dest): return dest # If already exists.
    url = 'http://pfam.sanger.ac.uk/family/%s/alignment/%s/gzipped' \
        % (pfam_family_id, atype)
    try:
        dest, hdrs = urllib.urlretrieve(url, dest)
    except IOError, e:
        print "Can't retrieve %r to %r: %s" % (url, dest_folder, e)
        return None
    return dest

def decompressGzipFile(gzip_file_name, destfolder, file_ext='ann'):
    '''
    Takes in GZIP file from PFAM or PDB and decompresses it to an
    .ann or specified file extension.
    '''
    try:
        with gzip.open(gzip_file_name, 'rb') as f:
            file_content = f.read()
    except:
        return None
    dest_file = os.path.join(destfolder, os.path.splitext \
                             (os.path.basename(gzip_file_name))[0] + '.' + \
                             file_ext)
    fout = open(dest_file, 'w')
    fout.write(file_content)
    fout.close()
    return dest_file

def readPfamFile(file_name):
    '''
    Takes in PFAM Stockholm file and reads it, returning
    a pfamFile object.
    '''         
    obj_out = pfamFile(file_name)
    obj_out.parse()
    return obj_out

def grabNCBIAccessionMetadata(accession_id):
    '''
    Takes in an accession ID and returns NCBI
    metadata from GenBank.
    '''        
    url = 'http://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?' + \
        'tool=portal&sendto=on&log$=seqview&db=protein&dopt=genpept' + \
        '&val=%s&extrafeat=984&maxplex=1' % (accession_id)
    return urllib.urlopen(url).read()

def grabPDBFile(pdb_code, dest_folder):
    '''
    Takes in a PDB code and acquires the PDB file
    from the PDB database, placing it in the specified
    destination folder.
    '''           
    dest = os.path.join(dest_folder, pdb_code + '.pdb')    
    if os.path.isfile(dest): return dest    
    url = 'http://www.rcsb.org/pdb/download/downloadFile.do?' + \
        'fileFormat=pdb&compression=NO&structureId=%s' \
        % (pdb_code)
    try:
        dest, hdrs = urllib.urlretrieve(url, dest)
    except IOError, e:
        print "Can't retrieve PDB for %s. Check network connection." % (pdb_code)
        return None        
    return dest

def printListToFile(list_in, dest):
    '''
    Prints a list to a file.
    '''           
    fout = open(dest, 'w')
    for item in list_in:
        fout.write(str(item) + '\n')
    fout.close()

def list2txttable(list_in, title):
    '''
    Converts a list to a string table.
    '''       
    txttable = texttable.Texttable()
    header = [title]
    txttable.header(header)
    txttable.set_deco(txttable.BORDER | txttable.HEADER | txttable.VLINES)
    txttable.set_chars(['-','|','+','-'])
    for item in list_in:
        txttable.add_row([item])
    return txttable.draw()

def removeListDuplicates(list_in): 
    '''
    Removes all duplicates in a list.
    '''     
    seen = {}
    list_out = []
    for item in list_in:
        if item in seen: continue
        seen[item] = 1
        list_out.append(item)
    return list_out
