# io.py
# -------------------------
# May 29, 2012; Alex Safatli
# -------------------------
# Helper functions that work
# with I/O operations and OS
# operations.

import os, glob
import shutil

def getFileName(path):
    '''
    Returns a filename of a given path.
    '''
    return os.path.splitext(os.path.basename(path))[0]

def getFolderName(path):
    '''
    Returns the directory name of a given path.
    '''
    return os.path.dirname(path)

def makeFolder(folder_name):
    '''
    Makes a folder. Returns path.
    '''
    if not os.path.exists(folder_name): 
        os.makedirs(folder_name)
    return folder_name

def moveFile(origin,target):
    '''
    Moves a file. Returns path.
    '''
    shutil.move(origin,target)
    return os.path.join(target,os.path.basename(origin))

def copyFile(origin,target):
    '''
    Copies a file. Returns path.
    '''
    shutil.copy(origin,target)
    return os.path.join(target,os.path.basename(origin))

def removeFile(origin):
    '''
    Removes a file.
    '''
    if os.path.isfile(origin):
        os.remove(origin)
        
def getFilesInFolder(folder, ext):
    '''
    Returns a list of all full file paths
    for a given extension in a given folder.
    '''
    return glob.glob(os.path.join(folder, '*.' + ext))  
