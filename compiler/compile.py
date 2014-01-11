#!/usr/bin/python

import sys
import os
import getopt
import shutil

# The order in which to compile the files
FILE_LIST = [

    'phys.js'

]

CLOSURE_COMPILER_FILENAME = 'compiler.jar'

VERBOSE_WARNINGS = True

OFF_OPTIONS = [
    'globalThis',
    'nonStandardJsDocs',
    'checkTypes'
]

WARNING_OPTIONS = [
    'accessControls',
    'checkDebuggerStatement',
    'checkRegExp',
    'const',
    'constantProperty',
    'strictModuleDepCheck',
    'visibility'
]

DEBUG = False                   # Whether or not to print verbose information

TEMP_DIR = 'temp'               # The name of the directory to place temporary file into

REMOVE_TEMP_DIR = True          # Whether or not to delete the temporary directory after compilation



def compileJS(filenames, outputFileName, compilationLevel, definitions=[]):

    if DEBUG:
        print 'Compiling code'
        print 'Input directory: '+INPUT_DIR
        print 'Output file: '+outputFileName
        print 'Temporary directory: '+TEMP_DIR
        print 'Tags to define: '+str(definitions)
        print 'Compilation Level: '+compilationLevel
        print

    # Preprocess

    cmd = 'python '+jsmacroFile+' --srcdir '+INPUT_DIR+' --dstdir '+TEMP_DIR
    if len(definitions):
        cmd += ' --def '+' --def '.join(definitions)

    os.system(cmd)

    # Compile the javascript

    if DEBUG:
        print 'Compiling into file "'+outputFileName+'"'

    compileCommand = 'java -jar '+CLOSURE_COMPILER_FILENAME+' --compilation_level '+compilationLevel

    if VERBOSE_WARNINGS:
        compileCommand += ' --warning_level=VERBOSE'

    for option in OFF_OPTIONS:
        compileCommand += ' --jscomp_off='+option

    for option in WARNING_OPTIONS:
        compileCommand += ' --jscomp_warning='+option

    filenames = os.listdir(TEMP_DIR)
    for filename in FILE_LIST:
        if filename in filenames:
            compileCommand += ' --js='+os.path.join(TEMP_DIR, filename)
            filenames.remove(filename)
        else:
            print('Did not find expected file "'+filename+'"')

    if len(filenames):
        print('Found extra files: '+','.join(filenames))

    compileCommand += ' --js_output_file='+outputFileName

    if DEBUG:
        print '    '+compileCommand
        print

    os.system(compileCommand)



if __name__ == '__main__':

    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "h:s:d:",
            ["help", "jsmacro=", "srcdir=", "dstdir="]
        )
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)

    # Handle parameters

    INPUT_DIR = '../source'         # The directory to compile the code from

    OUTPUT_DIR = '../'

    # Get the version

    version = 'latest'

    jsmacroFile = 'jsmacro.py'

    for o, a in opts:

        if o == "-h" or o == "--help":
            print("TODO")
            sys.exit(0)

        if o == "--jsmacro":
            jsmacroFile = a

        if o == "--srcdir":
            INPUT_DIR = a

        if o == "--dstdir":
            OUTPUT_DIR = a

        if o == "-version" or o == "-v":
            version = a

    # Create the temporary directory

    if DEBUG:
        print
        print 'Creating temporary directory "'+TEMP_DIR+'"'

    if os.path.isdir(TEMP_DIR):
        shutil.rmtree(TEMP_DIR)

    os.mkdir(TEMP_DIR)

    # Get the filenames

    if DEBUG:
        print 'Reading files in from directory "'+INPUT_DIR+'"'

    filenames = os.listdir(INPUT_DIR)

    # Compile multiple versions of the code
    compileJS(filenames, OUTPUT_DIR+'phys.min.js', 'SIMPLE_OPTIMIZATIONS')

    # Remove the temporary directory
    if REMOVE_TEMP_DIR:
        shutil.rmtree(TEMP_DIR)
