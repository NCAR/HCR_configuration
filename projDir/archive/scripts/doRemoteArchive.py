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
        print >>sys.stderr, "Running: ", os.path.basename(__file__)
        print >>sys.stderr, "  Options:"
        print >>sys.stderr, "    Debug: ", options.debug
        print >>sys.stderr, "    Verbose: ", options.verbose
        print >>sys.stderr, "    Remote host: ", options.host
        print >>sys.stderr, "    Source dir: ", options.dir
        print >>sys.stderr, "    Project name: ", options.projectName
        print >>sys.stderr, "    Drive index: ", driveIndex
        
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
        print >>sys.stderr, "======================="
        print >>sys.stderr, "Target disk drive list:"
        for drive in driveList:
            print >>sys.stderr, \
                  "  drive, device: ", \
                  drive, deviceTable[drive]
        print >>sys.stderr, "Drive to use: ", driveToUse

    if (driveToUse == "none"):
        print >>sys.stderr, "ERROR - no drive at index: ", driveIndex
        exit(1)

    # compute day string for today

    now = time.gmtime()
    nowTime = datetime.datetime(now.tm_year, now.tm_mon, now.tm_mday,
                                now.tm_hour, now.tm_min, now.tm_sec)
    nowStr = nowTime.strftime("%Y%m%d")
    nowSecs = time.mktime(nowTime.timetuple())

    # compute the earliest valid time

    if (options.debug == True):
        print >>sys.stderr, "======================="
        print >>sys.stderr, "Time details: "
        print >>sys.stderr, "  now time: ", nowTime
        print >>sys.stderr, "      secs: ", nowSecs

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
        tokens = line.split()
        if (tokens[0].find('/dev') >= 0):
            partition = tokens[5]
            if (partition.find('HCR') >= 0):
                driveList.append(partition)
                deviceName = tokens[0]
                deviceTable[partition] = deviceName

    # sort drive list
    
    driveList.sort()

########################################################################
# archive specified dir to specified drive

def doArchiveToDrive(drive):

    if (options.debug == True):
        print >>sys.stderr, "====================================="
        print >>sys.stderr, "Syncing host: ", options.host
        print >>sys.stderr, "         dir: ", options.dir
        print >>sys.stderr, "    to drive: ", drive

    # compute target dir
    
    targetDir = drive
    targetDir = os.path.join(targetDir, "hcr")
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
        print >>sys.stderr, "running cmd:",cmd
    
    try:
        retcode = subprocess.call(cmd, shell=True)
        if retcode < 0:
            print >>sys.stderr, "Child was terminated by signal: ", -retcode
        else:
            if (options.debug == True):
                print >>sys.stderr, "Child returned code: ", retcode
    except OSError, e:
        print >>sys.stderr, "Execution failed:", e

########################################################################
# Run - entry point

if __name__ == "__main__":
   main()
