from __future__ import print_function
import os
import sys
from subprocess import Popen, PIPE

def CheckFileName(files, badChars=['\\','/',':','*','?','\"','<','>','|'],
                  replacement=None):
    """
    Check if file does not have \ / : * ? \" < > |. If given a replacement
    string, replace the erroneous character with it, and return a list of fixed
    file names. If a replacement is not given, an Exception is thrown instead.
    
    Note this is for Python2, if files is unicode, then it will return a list.
    Of course in Python3, every string is 'unicode' while the unicode keyword
    doesn't exist anymore. So if you're not sure, just use:
    CheckFileName(str(files), ...) if you're using Python2 and 'files' is either
    an ASCII string or unicode string.
    """

    # If the 'files' is a string, remember that fact and make it a list
    isString = isinstance(files, str)
    if isString:
        files = [files]

    # Eliminate erroneous characters
    result = []
    for fileName in files:
        for char in badChars:
            pos = fileName.find(char)
            while pos != -1: # Erroneous character
                if replacement is None:
                    raise Exception("File name contains illegal character: %c"
                                    % char)
                fileName = fileName[:pos]+replacement+fileName[pos+len(char):]
                pos = fileName.find(char)
        result.append(fileName)

    # If a string was given, return a string instead of a list
    if isString:
        return result[0]
    else:
        return result

def CheckFiles(*files):
    """Check if files exists"""
    for file_name in files:
        if file_name == None:
            #I may find a way to use another parameter that isn't finicky
            #i.e. only Python3 supports def CheckFiles(*files, exceptions=False)
            #if exceptions:
            #  raise Exception("File is labelled as None, missing parameter?\n")
            #else:
            #    return False
            return False
        if not os.path.isfile(file_name):
            #if exceptions:
            #    raise Exception("File does not exist: %s\n" % file_name)
            #else:
            #    return False
            return False
    return True

def CheckModule(module):
    """Raises ImportError if module does not exist"""
    version = sys.version_info.major
    subv = sys.version_info.minor
    if version == 2:
        import imp
        imp.find_module(module)
    else: #version == 3
        import importlib
        if subv <= 3:
            if not importlib.find_loader(module):
                raise ImportError("No module named %s" % module)
        else: #subv > 4
            if not importlib.util.find_spec(module):
                raise ImportError("No module named %s" % module)

def CheckParam(param, test=(lambda x: x), stream=sys.stderr, fail_str=None,
               exit=False):
    """See if param passes test, else print fail_str and possibly exit"""
    if not test(param):
        if fail_str:
            FPrint(stream, fail_str)
        if exit:
            sys.exit(2)
        return False
    return True

def CleanUp(*files):
    """Deletes all files given"""
    for file_name in files:
        Command("rm", "-f", file_name)

def Command(*args):
    """Returns stdout, throws stderr"""
    stdout, stderr = Popen(list(args), stdout=PIPE, stderr=PIPE).communicate()
    if stderr != "":
        raise Exception(stderr)
    return stdout            
            
def CreateDir(dir_name, reuse=True):
    """Creates a directory and returns name of the directory"""
    if os.path.isdir(dir_name):
        if reuse == False:
            count = 1
            while os.path.isdir(dir_name + ("_%i" % count)):
                count += 1
            dir_name += ("_%i" % count)
            Command("mkdir", dir_name)
    else:
        Command("mkdir", dir_name)
    return dir_name

def FPrint(stream, *args, **kwargs):
    """Prints to stream (either sys.stdout, sys.stderr, or file)"""
    print(*args, file=stream, **kwargs)
    
