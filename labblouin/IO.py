''' 
Different shortcut functions for dealing with file systems.

IO Python Library / 2012 / Alex Safatli

E-mail: safatli@cs.dal.ca
Dependencies: -

'''

from os.path import splitext, basename, dirname, exists, join, isfile, isdir
from os import makedirs, remove
import glob
import shutil

def getFileName(path):

    ''' Returns a filename of a given path. '''

    return splitext(basename(path))[0]

def withExtension(file_name,ext,basename=True):

    ''' Add an extension to a file name. '''

    if not basename: file_name = getFileName(file_name)
    return file_name + "." + ext

def getFolderName(path):

    ''' Returns the directory name of a given path. '''

    return dirname(path)

def joinFileName(file_path,join_path):

    ''' Gets the filename of file_path and joins it, as os.path.join,
    to join_path. '''

    return join(getFileName(file_path),join_path)

def makeFolder(folder_name):

    ''' Makes a folder (if does not already exist). Returns path. '''

    if not exists(folder_name): makedirs(folder_name)
    return folder_name

def makeFolderFromFileName(file_path):

    ''' Makes a folder name equal to the base file name of the given file_path
    as os.path.basename. '''

    return makeFolder(getFileName(file_path))

def moveFile(origin,target):

    ''' Moves a file. Returns path. '''

    shutil.move(origin,target)
    return join(target,basename(origin))

def copyFile(origin,target):

    ''' Copies a file. Returns path. '''

    shutil.copy(origin,target)
    return join(target,basename(origin))

def removeFile(origin):

    ''' Removes a file. '''

    if isfile(origin): remove(origin)
        
def getFilesInFolder(folder, ext):

    ''' Returns a list of all full file paths
    for a given extension in a given folder. '''

    return glob.glob(join(folder, '*.' + ext))  
