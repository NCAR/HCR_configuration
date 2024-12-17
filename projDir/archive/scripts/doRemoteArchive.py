#!/usr/bin/env python

# ========================================================================== #
#
# Archive dir from remote host - most likely rds
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
    global driveList
    global deviceTable
    global nowSecs
    global driveIndex

    # set some variables from the environment

    dataDir = os.environ['DATA_DIR']
    projDir = os.environ['PROJ_DIR']
    projName = os.environ['PROJECT_NAME']
    
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
    parser.add_option('--host',
                      dest='host', default='rds',
                      help='Name of remote host')
    parser.add_option('--dir',
                      dest='dir', default='projDir/logs',
                      help='Path of source directory')
    parser.add_option('--project',
                      dest='projectName', default=projName,
                      help='Name of project - used to generate target dir')
    parser.add_option('--driveIndex',
                      dest='driveIndex', default='0',
                      help='Which drive to use [0 or 1]')
    
    (options, args) = parser.parse_args()

    driveIndex = int(options.driveIndex)

    if (options.verbose == True):
        options.debug = True

    if (options.debug == True):
        print("Running: ", os.path.basename(__file__), file=sys.stderr)
        print("  Options:", file=sys.stderr)
        print("    Debug: ", options.debug, file=sys.stderr)
        print("    Verbose: ", options.verbose, file=sys.stderr)
        print("    Remote host: ", options.host, file=sys.stderr)
        print("    Source dir: ", options.dir, file=sys.stderr)
        print("    Project name: ", options.projectName, file=sys.stderr)
        print("    Drive index: ", driveIndex, file=sys.stderr)
        
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

    # compute day string for today

    now = time.gmtime()
    nowTime = datetime.datetime(now.tm_year, now.tm_mon, now.tm_mday,
                                now.tm_hour, now.tm_min, now.tm_sec)
    nowStr = nowTime.strftime("%Y%m%d")
    nowSecs = time.mktime(nowTime.timetuple())

    # compute the earliest valid time

    if (options.debug == True):
        print("=======================", file=sys.stderr)
        print("Time details: ", file=sys.stderr)
        print("  now time: ", nowTime, file=sys.stderr)
        print("      secs: ", nowSecs, file=sys.stderr)

    # perform the archival, to selected drive
        
    doArchiveToDrive(driveToUse)

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
# archive specified dir to specified drive

def doArchiveToDrive(drive):

    if (options.debug == True):
        print("=====================================", file=sys.stderr)
        print("Syncing host: ", options.host, file=sys.stderr)
        print("         dir: ", options.dir, file=sys.stderr)
        print("    to drive: ", drive, file=sys.stderr)

    # compute target dir
    
    targetDir = drive
    targetDir = os.path.join(targetDir, "rsf")
    targetDir = os.path.join(targetDir, "archive")
    targetDir = os.path.join(targetDir, options.projectName)
    targetDir = os.path.join(targetDir, options.host)
    targetDir = os.path.join(targetDir, options.dir)
    
    # ensure target dir exists

    cmd = "mkdir -p "
    cmd += targetDir
    runCommand(cmd)

    # compute rsync command

    cmd = "rsync -av --rsh='ssh' "
    cmd += options.host
    cmd += ":"
    cmd += options.dir
    cmd += "/\* "
    cmd += targetDir
    # for now, run in foreground
    # cmd += " &"

    # run the command
    
    runCommand(cmd)

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
            if (options.debug == True):
                print("Child returned code: ", retcode, file=sys.stderr)
    except OSError as e:
        print("Execution failed:", e, file=sys.stderr)

########################################################################
# Run - entry point

if __name__ == "__main__":
   main()
