import Bio
from Bio.PDB import PDBList, PDBIO, PDBParser
from Bio.PDB.PDBIO import Select
from labblouin import homology
from io import StringIO
import pandas as pd
from multiprocessing import Pool
from subprocess import Popen, PIPE
import optparse
import csv
import os
#from PDBnet import PDBstructure as st


def iterfasta(fasta):
    """
    parse a fasta and yield the info
    :param fasta: fata file
    :return: iterator
    """
    with open(fasta) as F:
        spl = F.read().strip().split('>')
        for seq in spl:
            if seq == '':
                continue
            else:
                sp = seq.split('\n')
                yield sp[0], '\n'.join(sp[1:]).strip()


def uniblast(sequence, db, evalue, tgts):
    os.environ['BLASTDB'] = db
    seq ='>%s\n%s' % sequence
    blast_line = ['psiblast', '-db', db, '-evalue', str(evalue),
                  '-num_iterations', '0', '-max_target_seqs', str(tgts),
                  '-outfmt', '6 qaccver saccver pident evalue qcovs length '
                             'staxid']
    taxon_line = ['taxonkit', 'lineage', '-i', '7']
    taxon_line2 = ['taxonkit', 'reformat', '-i', '8', '-f', '{s}']
    blast = Popen(blast_line, stdin=PIPE, stdout=PIPE, stderr=PIPE,
                  env=os.environ.copy(), encoding='utf8', cwd=os.getcwdu())
    o, e = blast.communicate(seq)
    txnkit1 = Popen(taxon_line, stdin=PIPE, stdout=PIPE,  stderr=PIPE,
                    env=os.environ.copy(), encoding='utf8', cwd=os.getcwdu())
    o1, e1 = txnkit1.communicate(o[:o.find('\n\n')] )
    txnkit2 = Popen(taxon_line2, stdin=PIPE, stdout=PIPE, stderr=PIPE,
                    env=os.environ.copy(), encoding='utf8', cwd=os.getcwdu())
    o2, e2 = txnkit2.communicate(o1)
    df = pd.read_table(StringIO(o2), sep='\t', header=None,
                       names='qaccver saccver pident evalue qcovs length '
                             'staxid lineage species'.split(),
                       quoting=csv.QUOTE_NONE, encoding='utf-8')
    df = df.sort_values(by=['pident', 'evalue', 'qcovs', 'length'],
                        ascending=[False, True, False, False])
    df = df[(df.pident > 30) & (df.qcovs > 70)][~df.saccver.duplicated()]
    pdbl = PDBList()
    pdb_codes = df.saccver.unique()
    for i in pdb_codes:
        pdb = i[:4]
        ch = i[-1]
        # Download pdb
        file_path = pdbl.retrieve_pdb_file(pdb, pdir='PDBs', file_format='pdb')
        parser = PDBParser()
        st = parser.get_structure(pdb, file_path)
        ou = PDBIO()
        ou.set_structure(st)
        os.remove(os.path.join('PDBs', pdb))
    return df
