#!/usr/bin/env python

# ========================================================================== #
#
# Configure the host for a given role
#
# ========================================================================== #

from __future__ import print_function
import os
import sys
from optparse import OptionParser
import datetime
import subprocess
import string

def main():

    global options

    homeDir = os.environ['HOME']
    projDir = os.path.join(homeDir, 'projDir')
    controlDir = os.path.join(projDir, 'control')
    defaultGitDir = os.path.join(homeDir, "git/HCR_configuration")

    # parse the command line

    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option('--debug',
                      dest='debug', default=True,
                      action="store_true",
                      help='Set debugging on')
    parser.add_option('--verbose',
                      dest='verbose', default=False,
                      action="store_true",
                      help='Set verbose debugging on')
    parser.add_option('--gitDir',
                      dest='gitDir', default=defaultGitDir,
                      help='Path of main directory in git')
    parser.add_option('--dataDir',
                      dest='dataDir', default='/data/hcr',
                      help='Path of installed data dir')
    (options, args) = parser.parse_args()
    
    if (options.verbose):
        options.debug = True

    # compute paths

    gitProjDir = os.path.join(options.gitDir, 'projDir')
    gitSystemDir = os.path.join(gitProjDir, 'system')
    
    # debug print

    if (options.debug):
        print >>sys.stderr, "Running script: ", os.path.basename(__file__)
        print >>sys.stderr, ""
        print >>sys.stderr, "  Options:"
        print >>sys.stderr, "    Debug: ", options.debug
        print >>sys.stderr, "    Verbose: ", options.verbose
        print >>sys.stderr, "    homeDir: ", homeDir
        print >>sys.stderr, "    projDir: ", projDir
        print >>sys.stderr, "    controlDir: ", controlDir
        print >>sys.stderr, "    gitDir: ", options.gitDir
        print >>sys.stderr, "    gitProjDir: ", gitProjDir
        print >>sys.stderr, "    gitSystemDir: ", gitSystemDir

    # read current host type if previously set

    prevHostType = 'archiver'
    hostTypePath = os.path.join(homeDir, '.host_type')
    if (os.path.exists(hostTypePath)):
        hostTypeFile = open(hostTypePath, 'r')
        prevHostType = hostTypeFile.read()
        prevHostType = prevHostType.strip(string.whitespace)
    if (options.debug):
        print >>sys.stderr, "    prevHostType: ", prevHostType

    # get the host type interactively

    hostTypes = [ 'archiver', 'drx' ]

    print("", file=sys.stdout)
    print("Choose host type from the following list", file=sys.stdout)
    print("       or hit enter to use host type shown:", file=sys.stdout)
    for hostType in hostTypes:
        print("     ", hostType, file=sys.stdout)
    if (sys.version_info > (3, 0)):
        hostType = input('    ............. (' + prevHostType + ')? ')
    else:
        hostType = raw_input('    ............. (' + prevHostType + ')? ')
    if (len(hostType) < 4):
        hostType = prevHostType
    else:
        typeIsValid = False
        for htype in hostTypes:
            if (hostType == htype):
                typeIsValid = True
        if (typeIsValid != True):
            print("ERROR - invalid host type: ", hostType, file=sys.stderr)
            sys.exit(1)

    gitProjDir = os.path.join(options.gitDir, "projDir")

    # save the host type to ~/.host_type

    hostTypeFile = open(hostTypePath, "w")
    hostTypeFile.write(hostType + '\n')
    hostTypeFile.close()

    # banner

    print(" ", file=sys.stdout)
    print("*********************************************************************", file=sys.stdout)
    print()
    print("  configure HCR", file=sys.stdout)
    print()
    print("  runtime: " + str(datetime.datetime.now()), file=sys.stdout)
    print()
    print("  host type: ", hostType, file=sys.stdout)
    print()
    print("*********************************************************************", file=sys.stdout)
    print(" ", file=sys.stdout)

    # make links to the dotfiles in git projDir
    
    os.chdir(homeDir)
    for rootName in ['cshrc', 'bashrc', 'emacs', 'Xdefaults' ]:
        dotName = '.' + rootName
        removeSymlink(homeDir, dotName)
        sourceDir = os.path.join(gitSystemDir, 'dotfiles')
        sourcePath = os.path.join(sourceDir, rootName)
        cmd = "ln -s " + sourcePath + " " + dotName
        runCommand(cmd)

    # make link to projDir

    removeSymlink(homeDir, 'projDir')
    os.chdir(homeDir)
    cmd = "ln -s " + gitProjDir
    runCommand(cmd)
    
    # make link to proc_list, crontab and data_list

    os.chdir(controlDir)
    cmd = "ln -s proc_list." + hostType + " proc_list"
    runCommand(cmd)
    cmd = "ln -s crontab." + hostType + " crontab"
    runCommand(cmd)
    cmd = "ln -s data_list." + hostType + " data_list"
    runCommand(cmd)
    
    ############################################
    # data dir - specific to the host type
    # populate installed data dir /data/hcr
    
    dataDirsPath = os.path.join(options.gitDir, 'data_dirs')
    dataSubDir = "data." + hostType
    templateDataDir = os.path.join(dataDirsPath, dataSubDir)
    installDataDir = os.path.join(options.dataDir, dataSubDir)

    if (options.debug):
        print("Install data dir: ", installDataDir, file=sys.stderr)

    # create symlink to data if it does not already exist

    os.chdir(projDir)
    if (os.path.exists('data') == False):
        cmd = "ln -s " + installDataDir + " data"
        runCommand(cmd)

    # create symlink to logs

    os.chdir(projDir)
    if (os.path.exists('logs') == False):
        cmd = "ln -s data/logs"
        runCommand(cmd)

    # rync template dir into data dir

    os.chdir(templateDataDir)
    cmd = "rsync -av * " + installDataDir
    runCommand(cmd)

    # create symlinks for parameter files
    # from data tree back into the template

    debugStr = ""
    if (options.debug):
        debugStr = " --debug"
    cmd = "createParamLinks.py --templateDir " + \
          templateDataDir + " --installDir " + installDataDir + debugStr
    runCommand(cmd)
    
    # done

    sys.exit(0)
    
########################################################################
# Remove a symbolic link
# Exit on error

def removeSymlink(dir, linkName):

    os.chdir(dir)

    # remove if exists
    if (os.path.islink(linkName)):
        os.unlink(linkName)
        return

    if (os.path.exists(linkName)):
        # link name exists but is not a link
        print("ERROR - trying to remove symbolic link", file=sys.stderr)
        print("  dir: ", dir, file=sys.stderr)
        print("  linkName: ", linkName, file=sys.stderr)
        print("This is NOT A LINK", file=sys.stderr)
        sys.exit(1)

########################################################################
# Run a command in a shell, wait for it to complete

def runCommand(cmd):

    if (options.debug == True):
        print("running cmd:",cmd, file=sys.stderr)

    try:
        retcode = subprocess.call(cmd, shell=True)
        if retcode < 0:
            print("Child was terminated by signal: ", -retcode, file=sys.stderr)
        else:
            if (options.verbose == True):
                print("Child returned code: ", retcode, file=sys.stderr)
    except OSError as e:
        print("Execution failed:", e, file=sys.stderr)

########################################################################
# Run - entry point

if __name__ == "__main__":
   main()
