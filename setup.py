''' An installation script for these tools. '''

# Date:   Oct 17 2014
# Author: Alex Safatli
# E-mail: safatli@cs.dal.ca

from setuptools import setup

VERSION = '0.3.1'
DESCRIP = 'Bioinformatic python utilities modules, and libraries for use with Lab Blouin applications.'
URL     = 'http://www.github.com/LabBlouin/LabBlouinTools'
AUTHOR  = 'Christian Blouin et al.'
EMAIL   = 'cblouin@cs.dal.ca'
DEPNDS  = ['scikit-learn']
LINKS   = ['svn+https://pymol.svn.sourceforge.net/svnroot/pymol/trunk/pymol']

setup(name='labblouintools',version=VERSION,description=DESCRIP,url=URL,author=AUTHOR,author_email=EMAIL,license='GNU',packages=['labblouin'],install_requires=DEPNDS,dependency_links=LINKS,zip_safe=False)
