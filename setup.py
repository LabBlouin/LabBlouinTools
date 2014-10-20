''' An installation script for these tools. '''

# Date:   Oct 17 2014
# Author: Alex Safatli
# E-mail: safatli@cs.dal.ca

from setuptools import setup

VERSION = '0.1'
DESCRIP = 'Bioinformatic python utilities modules, and libraries for use with Lab Blouin applications.'
URL     = 'http://www.github.com/AlexSafatli/LabBlouinTools'
AUTHOR  = 'Alex Safatli'
EMAIL   = 'safatli@cs.dal.ca'
DEPNDS  = ['pymol']
LINKS   = []

setup(name='labblouintools',version=VERSION,description=DESCRIP,url=URL,author=AUTHOR,author_email=EMAIL,license='MIT',packages=['labblouin'],install_requires=DEPNDS,dependency_links=LINKS,zip_safe=False)
