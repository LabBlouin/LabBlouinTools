#!/bin/python

''' 
Different shortcut functions for dealing with file systems.

IO Python Library / 2012 / Alex Safatli

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

E-mail: safatli@cs.dal.ca
Dependencies: -

'''

import os, glob
import shutil

def getFileName(path):
    ''' Returns a filename of a given path. '''
    return os.path.splitext(os.path.basename(path))[0]

def getFolderName(path):
    ''' Returns the directory name of a given path. '''
    return os.path.dirname(path)

def makeFolder(folder_name):
    ''' Makes a folder (if does not already exist). Returns path. '''
    if not os.path.exists(folder_name): 
        os.makedirs(folder_name)
    return folder_name

def moveFile(origin,target):
    ''' Moves a file. Returns path. '''
    shutil.move(origin,target)
    return os.path.join(target,os.path.basename(origin))

def copyFile(origin,target):
    ''' Copies a file. Returns path. '''
    shutil.copy(origin,target)
    return os.path.join(target,os.path.basename(origin))

def removeFile(origin):
    ''' Removes a file. '''
    if os.path.isfile(origin):
        os.remove(origin)
        
def getFilesInFolder(folder, ext):
    ''' Returns a list of all full file paths
    for a given extension in a given folder. '''
    return glob.glob(os.path.join(folder, '*.' + ext))  
