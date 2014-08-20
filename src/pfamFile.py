#!/bin/python

# pfamFile.py
# -------------------------
# May 14, 2012; Alex Safatli
# -------------------------
# Handles all file manipulation
# and parsing for a PFAM flat
# file (Stockholm file format).

import math

class pfam:
    
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
        