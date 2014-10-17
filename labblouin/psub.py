#!/usr/bin/env python

''' Run as many of a input set of jobs/commands as there are cores at any given time (i.e. emulate a queue such as qsub present on GRID-powered hardware). '''

# Date:   Oct 2 2014
# Author: Alex Safatli
# E-mail: safatli@cs.dal.ca

import sys, os
from multiprocessing.pool import Pool

# A function to call on a worker.

def call(instr):
    os.system(instr)

# Get input.

if len(sys.argv) != 3:
    print 'Usage: %s file_with_instructions num_cores' % (sys.argv[0])
    exit(2)
fnst = sys.argv[1]
ncrs = int(sys.argv[2])

# Read all instructions.

cmds = []
o = open(fnst,'r')
for line in o: cmds.append(line.strip('\n'))
o.close()

# Run jobs.

if ncrs > os.cpu_count(): ncrs = os.cpu_count()
pool = Pool(ncrs)
results = [pool.apply_async(call,(cmd,)) for cmd in cmds]
pool.close()
for _ in results: _.get()