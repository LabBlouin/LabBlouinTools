#!/usr/bin/python
'''
Get the amino acid sequence from a landmark file, align this profile to more sequences if required
using muscle, and mask the membership vector into this alignment.
'''
#importing bit####################################################################################################
import sys,os, glob
# End importing###################################################################################################

#Some definitions##################################################################################################
def InferSingleLettterCode(threlettercode):
	'''
	Convert the amino acid three letter code, into a single letter code
	'''
	from tl.rename.case import transform_sentence_case
	aa = {'Gly':'G' , 'Ala':'A' , 'Val':'V' , 'Leu':'L' , 'Ile':'I' , 'Met':'M' , 'Phe':'F' ,
	      'Trp':'W' , 'Pro':'P' , 'Ser':'S' , 'Thr':'T' , 'Cys':'C' , 'Tyr':'Y' , 'Asn':'N' , 
	      'Gln':'Q' , 'Asp':'D' , 'Glu':'E' , 'Lys':'K' , 'Arg':'R' , 'His':'H'}
	singleletter = aa[transform_sentence_case([threlettercode])[0]]
	return singleletter
	
	
def InferSeqs_landmarks(prefix):
	'''
	Given a landmark file, return a dictionary with the name as key and the
	sequence as value
	'''
	seqs={}
	inf = open(prefix+'.landmarks')
	rf = inf.read()
	spl = rf.split('\n>')
	for i in range(len(spl)):
		s=''
		name=spl[i][:spl[i].find('\n')]
		bline = spl[i][spl[i].find('\n')+1:].split()

		for j in range(2,len(bline),3):
			s+= InferSingleLettterCode(bline[j])
		seqs[name]=s
	return seqs


def dict2fasta(outfilename, dic):
	'''
	Convert the dictionary return by InferSeqs_landmarks fuction
	and writes a fasta file. The dictionary should have the seq name as
	key and a string as value containing the sequence
	'''
	fout=open(outfilename,'w')
	for k,v in dic.iteritems():
		fout.write('>'+k+'\n'+v+'\n')
	fout.close()

def fasta2dict(fastafilename):
	'''
	convert a fasta file into a python dictionary, being the keys the name of the sequence
	and the values the sequence itself, in a list
	'''
	dic={}
	inf = open(fastafilename)
	rf = inf.read().split('\n>')
	for i in range(len(rf)):
		s=''
		name=rf[i][rf[i].find('\n')-1]
		seq = rf[i][rf[i].find('\n'):].strip()
		dic[name]=str2list(seq)
	return dic
			
def Mem2List2tuple(prefix):
	'''
	Read the membership vector, and returns a list and a dictionary with tupples for
	each module, indicating the indexes
	'''
	mem = open(prefix+'.graphcluster').read()
	mem = mem.strip().split()
	s = set(mem)
	tuples = {}
	for i in range(len(s)):
		tuples[list(s)[i]]=()
		for j in range(len(mem)):
			if mem[j] == list(s)[i]:
				tuples[list(s)[i]]+=(j,)
		
	return mem, tuples

def str2list(string):
	'''
	returns a list of a string that have no spaces between elements
	'''
	l=[]
	for e in string:
		l.append(e)
	return l
	
def MaskAln(fastafilename, modulename, mem, tuples):
	'''
	Strip the module in an alignment (fasta) file
	'''
	ms=[]
	names=[]
	afadic = fasta2dict(fastafilename)
	for k, v in afadic.iteritems():
		s=''
		for e in tuples[modulename]:
			s+=v[e]
		names.append(k)
		ms.append(s)
	return ms, names
	
def modseq2Fasta(prefix,modseq,names):
	for k,v in modseq.iteritems():
		fout=open(prefix+'_%s.fst'%(k),'w')
		for e in range(len(v)):
			fout.write('>'+names[e]+'\n'+v[e]+'\n')
		fout.close()
	
# End of definitions###############################################################################################

# Aplication of the code ##########################################################################################
if __name__ == "__main__":
	# Some defaults
	seqname = None
	phylip = False
	# Get commands line
	for arg in sys.argv[1:]:
		if arg == '-align=':
			seqname = arg[7:]
		elif arg == '-phylip':
			phylip=True
	
	prefix = sys.argv[1]
	
	#create the fasta file of the structural alignment from which modules are being analized
	dict2fasta(prefix+'.fas', InferSeqs_landmarks(prefix))
	# if more sequences perform a profile analysis
	if seqname:
		os.system('muscle -profile -in1 %.fas -in2 %s.fa -out %s_combined.afa'%(prefix , seqname , prefix))
	# Otherwise, rename the fas file and use it as the aln.
	else:
		os.system('mv %s.fas %s.afa'%(prefix,prefix))
	
	# given the alignment and the membership vector, mask each module
	mem , tuples = Mem2List2tuple(prefix)
	modseq={}
	for e in list(set(mem)):
		if e == '?':
			continue
		else:
			modseq[e], names = MaskAln(prefix+'.afa', e, mem, tuples)
	
	modseq2Fasta(prefix,modseq,names)
	
	if phylip:
		cwd = os.getcwd()
		files = glob.glob(cwd+'/*.fst')
		for f in files:
			os.system('Fasta2Phylip %s %s.phy'%(f,f[:-4]))