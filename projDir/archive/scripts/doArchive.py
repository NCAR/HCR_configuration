#!/usr/bin/env python

# ========================================================================== #
#
# Trigger the archive script
#
# ========================================================================== #

import os
import sys
from optparse import OptionParser
import time
import datetime
from datetime import date
from datetime import timedelta
import subprocess

def main():

    global options
    global dirList
    global driveList
    global deviceTable
    global nowSecs
    global validStartSecs
    global maxAgeSecs
    global driveIndex

    # set some variables from the environment

    dataDir = os.environ['DATA_DIR']
    projDir = os.environ['PROJ_DIR']
    projName = os.environ['PROJECT_NAME']
    defaultDirListPath = \
        os.path.join(projDir, "archive/params/archiveDirList")
    
    # parse the command line

    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option('--debug',
                      dest='debug', default='False',
                      action="store_true",
                      help='Set debugging on')
    parser.add_option('--verbose',
                      dest='verbose', default='False',
                      action="store_true",
                      help='Set verbose debugging on')
    parser.add_option('--test',
                      dest='testMode', default='False',
                      action="store_true",
                      help='Use preset dates instead of today')
    parser.add_option('--source',
                      dest='sourceDir', default=dataDir,
                      help='Path of source directory')
    parser.add_option('--project',
                      dest='projectName', default=projName,
                      help='Name of project - used to generate target dir')
    parser.add_option('--dirListPath',
                      dest='dirListPath',
                      default=defaultDirListPath,
                      help='Path to file containing directory list')
    parser.add_option('--maxAgeHours',
                      dest='maxAgeHours', default='12',
                      help='Max age of files to be archived (hours)')
    parser.add_option('--driveIndex',
                      dest='driveIndex', default='0',
                      help='Which drive to use [0 or 1]')
    
    (options, args) = parser.parse_args()

    maxAgeSecs = int(float(options.maxAgeHours) * 3600.0)
    driveIndex = int(options.driveIndex)

    if (options.verbose == True):
        options.debug = True

    if (options.debug == True):
        print("Running: ", os.path.basename(__file__), file=sys.stderr)
        print("  Options:", file=sys.stderr)
        print("    Debug: ", options.debug, file=sys.stderr)
        print("    Verbose: ", options.verbose, file=sys.stderr)
        print("    Test mode: ", options.testMode, file=sys.stderr)
        print("    Source dir: ", options.sourceDir, file=sys.stderr)
        print("    Project name: ", options.projectName, file=sys.stderr)
        print("    Drive index: ", driveIndex, file=sys.stderr)
        print("    Dir list path: ", options.dirListPath, file=sys.stderr)
        print("    Max age (hours): ", options.maxAgeHours, file=sys.stderr)
        print("             (secs): ", maxAgeSecs, file=sys.stderr)
        
    # compile the list of target drives

    compileDriveList()

    # find the drive to use

    driveToUse = "none"
    index = 0
    for drive in driveList:
        if (driveIndex == index):
            driveToUse = drive
            break
        index = index + 1

    if (options.debug == True):
        print("=======================", file=sys.stderr)
        print("Target disk drive list:", file=sys.stderr)
        for drive in driveList:
            print("  drive, device: ", \
                  drive, deviceTable[drive], file=sys.stderr)
        print("Drive to use: ", driveToUse, file=sys.stderr)

    if (driveToUse == "none"):
        print("ERROR - no drive at index: ", driveIndex, file=sys.stderr)
        exit(1)

    # read in the directory list

    readDirList()

    if (options.debug == True):
        print("=======================", file=sys.stderr)
        print("Dir list:", file=sys.stderr)
        for dir in dirList:
            print("  -->> ", dir, file=sys.stderr)

    # compute day string for today
    # in test mode we set the daystr directly

    if (options.testMode == True):
        now = time.strptime("20110425010000", "%Y%m%d%H%M%S")
    else:
        now = time.gmtime()

    nowTime = datetime.datetime(now.tm_year, now.tm_mon, now.tm_mday,
                                now.tm_hour, now.tm_min, now.tm_sec)
    nowStr = nowTime.strftime("%Y%m%d")
    nowSecs = time.mktime(nowTime.timetuple())

    # compute the earliest valid time

    validPeriodDays = timedelta(float(options.maxAgeHours) / 24.0)
    validStartTime = nowTime - validPeriodDays
    validStartSecs = nowSecs - maxAgeSecs

    if (options.debug == True):
        print("=======================", file=sys.stderr)
        print("Time details: ", file=sys.stderr)
        print("  now time: ", nowTime, file=sys.stderr)
        print("      secs: ", nowSecs, file=sys.stderr)
        print("  files valid from: ", validStartTime, file=sys.stderr)
        print("              secs: ", validStartSecs, file=sys.stderr)

    # perform the archival, to selected drive
        
    for dir in dirList:
        doArchiveToDrive(driveToUse, dir)

    sys.exit(0)

########################################################################
# Get the list of drives

def compileDriveList():

    global driveList
    global deviceTable

    deviceTable = {}

    # run df to get drive stats
    
    pipe = subprocess.Popen('df', shell=True,
                            stdout=subprocess.PIPE).stdout
    lines = pipe.readlines()

    # load up drive list and associated tables
    
    driveList = []
    for line in lines:
        tokens = line.decode("utf-8").split()
        if (tokens[0].find('/dev') >= 0):
            partition = tokens[5]
            if (partition.find('RSF') >= 0):
                driveList.append(partition)
                deviceName = tokens[0]
                deviceTable[partition] = deviceName

    # sort drive list
    
    driveList.sort()

########################################################################
# Read in the directory list

def readDirList():

    global dirList

    fp = open(options.dirListPath, 'r')
    lines = fp.readlines()
    dirList = []
    for line in lines:
        if (line[0] != '#' and len(line) > 1):
            dirList.append(line.strip())


########################################################################
# archive specified dir to specified drive

def doArchiveToDrive(drive, dir):

    # compute source dir

    sourceDir = options.sourceDir
    sourceDir = os.path.join(sourceDir, dir)

    if (os.path.isdir(sourceDir) == False):
        print("=====================================", file=sys.stderr)
        print("Dir does not exist: ", sourceDir, file=sys.stderr)
        print("  skipping ...", file=sys.stderr)
        return

    if (options.debug == True):
        print("=====================================", file=sys.stderr)
        print("Syncing, dir: ", dir, file=sys.stderr)
        print("    to drive: ", drive, file=sys.stderr)
        print("   sourceDir: ", sourceDir, file=sys.stderr)

    # compute target dir
    
    targetDir = drive
    targetDir = os.path.join(targetDir, "rsf")
    targetDir = os.path.join(targetDir, "archive")
    targetDir = os.path.join(targetDir, options.projectName)
    targetDir = os.path.join(targetDir, dir)

    # compile the list of files to be archived

    processDir(sourceDir, targetDir, 0)

########################################################################
# process this dir

def processDir(sourceDir, targetDir, level):

    if (options.verbose == True):
        print("processDir: sourceDir, targetDir, level: ", \
            sourceDir, ", ", targetDir, ", ", level, file=sys.stderr)
    
    if (level > 5):
        # too deep
        print("ERROR - recursion too deep", file=sys.stderr)
        return

    # read through directory looking for files to archive

    fileList = []
    contents = os.listdir(sourceDir)

    for entry in contents:

        if (entry[0] == '_'):
            # ignore files starting with '_'
            continue

        entryPath = os.path.join(sourceDir, entry)

        if (os.path.isdir(entryPath)):
            
            # dir, so recurse

            newSourceDir = entryPath
            newTargetDir = os.path.join(targetDir, entry)

            if (options.verbose == True):
                print("  ===>> recursing to dir: ", newSourceDir, file=sys.stderr)

            processDir(newSourceDir, newTargetDir, level + 1)

        elif (os.path.isfile(entryPath)):

            # file, so add to list

            entryStat = os.stat(entryPath)
            fileModSecs = entryStat.st_mtime
            timeDiffSecs = nowSecs - fileModSecs

            if (options.verbose == True):
                print("  ===>> found path, modSecs, timeDiffSecs: ", \
                    entryPath, ", ", fileModSecs, ", ", timeDiffSecs, file=sys.stderr)

            # append to fileList

            if (fileModSecs > validStartSecs):
                if (options.verbose == True):
                    print("  ===>> appending file: ", entry, file=sys.stderr)
                fileList.append(entry)

    if (len(fileList) < 1):
        return

    # sort the list alphabetically
    
    #fileList.sort(lambda x, y: cmp(x.lower(),y.lower()))
    fileList.sort()
    
    if (options.debug == True):
        print("Full file list:", file=sys.stderr)
        for fileName in fileList:
            print("  ", fileName, file=sys.stderr)
        print("Only the most recent items will be rsync'd", file=sys.stderr)

    # create file list as string

    fileListStr = ""
    for fileName in fileList:
        fileListStr = fileListStr + fileName
        fileListStr = fileListStr + " "

    # ensure target dir exists

    cmd = "mkdir -p "
    cmd += targetDir
    runCommand(cmd)

    # compute base command

    cmd = "cd "
    cmd += sourceDir
    # NOTE --size-only is important when rsync-ing to FAT filesystems
    #cmd += "; rsync -av -6 --size-only"
    cmd += "; rsync -av -6 "
    cmd += fileListStr
    cmd += " "
    cmd += targetDir
    cmd += "; sync;"
    # for now, run in foreground
    # cmd += " &"

    # run the command
    
    runCommand(cmd)

    print("============ DONE ============", file=sys.stderr)
    

########################################################################
# compile file list archive

def compileFileList(topdir, subdir, level):

    global fileList
    if (level == 0):
        fileList = []

    if (options.verbose == True):
        print("compileFileList: topdir, subdir, level: ", \
            topdir, ", ", subdir, ", ", level, file=sys.stderr)
    
    if (level > 5):
        # too deep
        print("ERROR - recursion too deep", file=sys.stderr)
        return

    if (len(subdir) > 0):
        thisDir = os.path.join(topdir, subdir)
    else:
        # top level
        thisDir = topdir
        
    contents = os.listdir(thisDir)

    for entry in contents:

        if (entry[0] == '_'):
            # ignore files starting with '_'
            continue

        entryPath = os.path.join(thisDir, entry)
        entryRelPath = os.path.join(subdir, entry)

        if (os.path.isdir(entryPath)):
            
            # recurse
            if (options.verbose == True):
                print("  ===>> recursing to dir: ", entryPath, file=sys.stderr)

            compileFileList(topdir, entryRelPath, level + 1)

        elif (os.path.isfile(entryPath)):

            entryStat = os.stat(entryPath)
            fileModSecs = entryStat.st_mtime
            timeDiffSecs = nowSecs - fileModSecs

            if (options.verbose == True):
                print("  ===>> path, modSecs, timeDiffSecs: ", \
                    entryPath, ", ", fileModSecs, ", ", timeDiffSecs, file=sys.stderr)

            # append to fileList

            if (fileModSecs > validStartSecs):
                if (options.verbose == True):
                    print("  ===>> appending file: ", entryRelPath, file=sys.stderr)
                fileList.append(entryRelPath)

    return

########################################################################
# Run a command in a shell, wait for it to complete

def runCommand(cmd):

    if (options.verbose == True):
        print("running cmd:",cmd, file=sys.stderr)
    
    try:
        retcode = subprocess.call(cmd, shell=True)
        if retcode < 0:
            print("Child was terminated by signal: ", -retcode, file=sys.stderr)
        else:
            if (options.debug == True):
                print("Child returned code: ", retcode, file=sys.stderr)
    except OSError as e:
        print("Execution failed:", e, file=sys.stderr)

########################################################################
# Run - entry point

if __name__ == "__main__":
   main()
