# logfile.py
# --------------------------------
# logfile> May 6, 2013; Alex Safatli
# XMLfile> Jun 26, 2013; Alex Safatli
# --------------------------------
# Encapsulates logfile management
# as well as XML file reading and
# parsing for recording of discrete
# information.

import os, timer, re
from sys import stdout
from datetime import datetime
import xml.etree.cElementTree as xml

def parseTime(time):
    y,w,d,h,m,s = 0,0,0,0,0,0
    nums = re.split('[a-z]',time)[:-1]
    lers = re.split('[0-9]*',time)[1:]
    if 'y' in lers: y = int(nums[lers.index('y')])
    if 'w' in lers: w = int(nums[lers.index('w')])
    if 'd' in lers: d = int(nums[lers.index('d')])
    if 'h' in lers: h = int(nums[lers.index('h')])
    if 'm' in lers: m = int(nums[lers.index('m')])
    if 's' in lers: s = int(nums[lers.index('s')])
    return y, w, d, h, m, s

def calculateTime(timevector):
    y,w,d,h,m,s = timevector
    return y*3.15569e7+w*604800+d*86400+h*3600+m*60+s

def timeToString():
    return datetime.now().strftime("%y-%m-%d %H:%M")

class XMLfile():
    def __init__(self,xmlf,root='root'):
        self.path = xmlf
        self.root = xml.Element(root)
        self.tree = None
    def compile(self):
        self.tree = xml.ElementTree(self.root)
    def read(self):
        self.tree = xml.parse(self.path)
        self.root = self.tree.getroot()
    def write(self):
        self.compile()
        self.tree.write(self.path)
    def _add(self,ele,it,*args):
        out = xml.SubElement(ele,it)
        for arg in args:
            if len(arg) != 2: continue
            else:
                name,val = arg
                name = str(name)
                val  = str(val)
            if name == 'xml': out.text = val
            else: out.set(name,val)
        return out
    def add(self,ele,it,*args):
        out = self._add(ele,it,*args)
        self.write()
        return out

class logfile():

    def __init__(self, logf, enable=True, silent=False, overwrite=False):
        self.logfile = logf
        self.numat = 0
        self.totalnum = 100
        self.stime = timer.getTime()
        self.lastnum = self.numat
        self.enabled = enable
        self.silent  = silent
        if enable and (not os.path.isfile(logf) or overwrite):
            fout = open(logf,'w')
            fout.close()

    def setTotalNum(self,tn):
        self.totalnum = tn

    def incrementTimer(self, amt=1):
        self.numat += amt

    def updateTimer(self,numat):
        self.numat = numat

    # Format seconds.
    def __fs__(self,tot):
        y,w,d,h,m,s = 0,0,0,0,0,0
        if (tot < 1 or self.numat == 0): return '...'
        y = tot/(60*60*24*360)
        tot = tot % (60*60*24*360)
        w = tot/(60*60*24*7)
        tot = tot % (60*60*24*7)
        d = tot/(60*60*24)
        tot = tot % (60*60*24)
        h = tot/(60*60)
        tot = tot % (60*60)
        m = tot/(60)
        tot = tot % 60
        s = tot
        out = ''
        if (y >= 1): out += '%dy' % y
        if (w >= 1): out += '%dw' % w
        if (d >= 1): out += '%dd' % d
        if (h >= 1): out += '%dh' % h
        if (m >= 1): out += '%dm' % m
        if (s >= 1): out += '%ds' % s
        return out

    def _format(self,msg,status):
        prefix = '\033[1m'        
        if status == 'ERROR':
            prefix += '\033[91mERROR '
        elif status == 'NOTE':
            prefix += '\033[94mNOTE '
        elif status == 'WARNING':
            prefix += '\033[93mWARNING '
        return prefix + '\033[0m' + msg

    def _write(self,msg,temp=False,silent=None,progress=True,status=None):
        if silent == None: silent = self.silent
        if silent and not self.enabled: return
        numat  = self.numat+1
        totnum = self.totalnum+1
        tim    = timer.estimateRemainingTime(self.stime,numat,totnum) # time remaining
        cur    = timeToString()
        perc   = (float(numat)/totnum)*100
        suff   = '\r' if temp else ''
        fiout  = '!!%s!! %s' % (status,msg) if status else msg
        tmout  = self._format(msg, status) if status else msg
        if progress:
            fiout = '[%s | %3d%% | rem %10s] %s' % (cur,perc,self.__fs__(tim), fiout)
            tmout = '%s %3d%% rem %10s | %s' % (datetime.now().strftime("%H:%M"),perc,self.__fs__(tim), tmout)
        if not silent:
            stdout.write('\033[K' + tmout + suff)
            stdout.flush()
        if not self.enabled: return
        else:
            fout = open(self.logfile,'a')
            fout.write(fiout.strip('\n') + '\n')
            fout.close()

    def writeTemporary(self,msg,silent=None,status=None):
        self._write(msg,temp=True,silent=silent,status=status)

    def write(self,msg,silent=None,progress=True,status=None):
        self._write(msg + '\n',silent=silent,progress=progress,status=status)

    def writeElapsedTime(self):
        t = timer.getTime()
        self.write('Elapsed: %s' % (self.__fs__(t-self.stime)))
