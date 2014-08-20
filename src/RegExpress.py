#!/usr/bin/python
'''
Prosite pattern generator. Will turn an alignment into a regular expression
'''

#importing bit####################################################################################################
import sys
# End importing###################################################################################################

#some constants####################################################################################################
aa = {'A':'Ala','C':'Cys','D':'Asp','E':'Glu','F':'Phe','G':'Gly','H':'His','I':'Ile','K':'Lys','L':'Leu','M':'Met',\
      'N':'Asn','P':'Pro','Q':'Gln','R':'Arg','S':'Ser','T':'Thr','V':'Val','W':'Trp','Y':'Tyr'}
allaa=set('ACDEFGHIKLMNPQRSTVWY')
#Some definitions##################################################################################################

def check_if_is_align(dictionary_fasta):
	'''
	will check the read fasta to see if is an alignment
	'''
	lenght=[]
	for s in dictionary_fasta.itervalues():
		lenght.append(len(s))
		
	if len(set(lenght)) != 1:
		print 'Sequences not of the same lenght. An alignment should be provided.'
		sys.exit(-1)
	else:
		l=len(s)
		return l
def read_fasta(prefix):
	'''
	Wll read the fasta, place sequences in a dictionary. Will break if not an alignment
	'''
	fastas={}
	f = open(prefix+'.fasta').read().split('\n>')
	for e in f:
		if e == '':
			continue
		else:
			name=e[:e.find('\n')]
			seq=e[e.find('\n')+1:].strip().replace('\n','')
			fastas[name]=seq
	l= check_if_is_align(fastas)
	return fastas, l

def strip_gaps(matrix):
	'''
	given an alignment matrix, strip out the gap columns
	'''
	gapc=[]
	for e in matrix:
		if '-' in e:
			gapc.append(matrix.pop(matrix.index(e)))
	return gapc

def reg_express_col(col):
	'''
	given a column of an alignment, return the string corresponding to the regular expression of such column
	'''
	s = set(col)
	if s == allaa:
		re='x'
	elif len(s) >= 10 and s != allaa:
		ts=allaa.difference(s)
		re='{'
		for e in ts:
			re+=e
		re+='}'
	elif 1 < len(s) <= 10:
		re='['
		for e in s:
			re+=e
		re+=']'
	elif len(s) == 1:
		for e in s:
			re=e
	return re
		
		
def reg_express(fastas,l):
	'''
	will create a matrix with the alignment, discard columns with gaps, and create a regular expression similar to prosite
	'''
	m=[]
	for i in range(l):
		m.append([])
	seqn=0
	for v in fastas.itervalues():
		for e in range(len(v)):
			m[e].append(v[e])
	#gapc = strip_gaps(m)
	
	RE=''
	for col in m:
		re = reg_express_col(col)
		if '-' in re:
			continue
		else:
			RE+=re+'-'
	RE = RE[:-1]
	
	return RE	
			
# End of definitions###############################################################################################

# Aplication of the code ##########################################################################################		
if __name__ == "__main__":
	if len(sys.argv) == 1 or '-help' or '-h' in sys.argv:
		print 'Usage: RegExpress.py prefix'
	
	prefix = sys.argv[1]
	fastas , l = read_fasta(prefix)
	RE = reg_express( fastas , l )
	print RE