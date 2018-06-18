'''
GM2PDB
'''

import sys
from labblouin.PDBnet import PDBstructure

def parsePDB(prefix):
	D={}
	if isinstance(prefix, PDBstructure):
		if prefix.ismodel: things = prefix.models
		else: things = prefix.chains
		for t,l in things.iteritems():
			D[t] = str(l).split('\n')[1:]
	else:
		with open(prefix.replace('_gps','')+'.pdb') as F:
			for line in F:
				if line.startswith('MODEL'):
					name= line.strip().split()
					name =  name[1]
					D[name]=[]
				elif line.startswith('ATOM'):
					D[name].append(line)
	return D

def write_coor(lines,prefix):
	coor = open(prefix+'.coordinates','w')
	info = open(prefix+'.atominfo','w')
	for l in lines[0]:
		if l.startswith('TER'):
			continue
		else:
			rest = l[:27]
			info.write(rest+'\n')
	lastline = lines[0][-3]
	for k in sorted(lines.keys()):
		line = '>'+str(k)+';'
		#rest = ''
		for l in lines[k]:
			if 'ATOM' in l:
				bl = l.strip().split()
				x = bl[5]
				y = bl[6]
				z = bl[7]
				line += x+';'+y+';'+z+';'	
		coor.write(line+'\n')
	coor.close()
	info.close()
	return lastline

def AtomInfo(prefix):
	lines = []
	with open(prefix + '.atominfo','r') as F:
		for line in F:
			line = line.strip().split('\n')
			lines.append(line[0])
	return lines[:-2]

def Coord(GMfile):
	d = {}
	with open(GMfile,'r') as F:
		for line in F:
			line = line.strip().split(';')
			if line[-1] == '':
				del line[-1]
			try: modelname = float(line[0].strip('>'))
			except: 
				modelname = line[0].strip('"').strip('>model').strip('>Model').strip('>MODEL').strip('>')
				modelname = modelname[:modelname.find(':')]
				modelname = float(modelname)
			d[int(modelname)] = line[1:]
	return d
			
def WritePDB(output,lines,d,lastline):
	with open(output+'.pdb','wb') as fout:
		for k,v in d.iteritems():
			fout.write('MODEL        %s\n'%(k))
			num = 0
			for x in range(0,len(v)-1,3):				
				fout.write('%s    %7.3f %7.3f %7.3f  0.00  0.00\n'%(lines[num],float(v[x]),float(v[x+1]),float(v[x+2])))
				num += 1
			fout.write('%s\nENDMDL\n'%(lastline))

if __name__ == "__main__":	
	prefix = sys.argv[1]
	GMfile = sys.argv[2]
	output = sys.argv[3]
	
	pdb = parsePDB(prefix)
	lastline = write_coor(pdb, prefix)	
	lines = AtomInfo(prefix)
	d = Coord(GMfile)
	WritePDB(output, lines, d, lastline)
