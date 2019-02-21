#!/usr/bin/env python

#===========================================================================
#
# Produce plots for CMIGITS data for ARISTO
#
#===========================================================================

import os
import sys
import subprocess
from optparse import OptionParser
import numpy as np
from numpy import convolve
from numpy import linalg, array, ones
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import dates
import math
import datetime
import contextlib

def main():

    # globals

    global options
    global debug
    global startTime
    global endTime
    global figNum
    figNum = 0

    appName = 'PlotCmigits.aristo.py'
    projDir = os.environ['PROJ_DIR']
    dataDir = os.environ['DATA_DIR']
    global dataFilePath
    cmigitsFilePath = '/scr/snow2/rsfdata/projects/aristo-17/hcr/ascii/ac_georef/C-MIGITS_RF02.txt';
    sdn500FilePath = '/scr/snow2/rsfdata/projects/aristo-17/hcr/ascii/ac_georef/ARISTO.rf02.20170224_193000.to.20170224_223500.txt';

    # parse the command line

    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option('--debug',
                      dest='debug', default=False,
                      action="store_true",
                      help='Set debugging on')
    parser.add_option('--cmigits_file',
                      dest='cmigitsFilePath',
                      default=cmigitsFilePath,
                      help='File path for cmigits data')
    parser.add_option('--sdn500_file',
                      dest='sdn500FilePath',
                      default=sdn500FilePath,
                      help='File path for sdn500 data')
    parser.add_option('--widthMain',
                      dest='mainWidthMm',
                      default=300,
                      help='Width of main figure in mm')
    parser.add_option('--heightMain',
                      dest='mainHeightMm',
                      default=200,
                      help='Height of main figure in mm')
    parser.add_option('--filtLen',
                      dest='filtLen',
                      default=1,
                      help='Len of moving mean filter')
    parser.add_option('--decimation',
                      dest='decimation',
                      default=1,
                      help='Decimation factor when reading in data')
    parser.add_option('--start',
                      dest='startTime',
                      default='2017 02 24 19 30 00',
                      help='Start time for XY plot')
    parser.add_option('--end',
                      dest='endTime',
                      default='2017 02 24 22 30 00',
                      help='End time for XY plot')
    
    (options, args) = parser.parse_args()
    
    # time limits

    year, month, day, hour, minute, sec = options.startTime.split()
    startTime = datetime.datetime(int(year), int(month), int(day),
                                  int(hour), int(minute), int(sec))
    year, month, day, hour, minute, sec = options.endTime.split()
    endTime = datetime.datetime(int(year), int(month), int(day),
                                int(hour), int(minute), int(sec))

    if (options.debug == True):
        print >>sys.stderr, "  cmigitsFilePath: ", options.cmigitsFilePath
        print >>sys.stderr, "  sdn500FilePath: ", options.sdn500FilePath
        print >>sys.stderr, "  startTime: ", startTime
        print >>sys.stderr, "  endTime: ", endTime

    # read in column headers for cmigits data

    iret, cmigitsHdrs, cmigitsData = readCmigitsColumnHeaders(options.cmigitsFilePath)
    if (iret != 0):
        sys.exit(-1)

    # read in data for cmigits results

    cmigitsData, cmigitsTimes = readCmigitsInputData(options.cmigitsFilePath, cmigitsHdrs, cmigitsData)

    # load up the data arrays

    loadCmigitsDataArrays(cmigitsData, cmigitsTimes)

    # read in column headers for sdn500 data

    iret, sdn500Hdrs, sdn500Data = readSdn500ColumnHeaders(options.sdn500FilePath)
    if (iret != 0):
        sys.exit(-1)

    # read in data for sdn500 results

    sdn500Data, sdn500Times = readSdn500InputData(options.sdn500FilePath, sdn500Hdrs, sdn500Data)

    # load up the data arrays

    loadSdn500DataArrays(sdn500Data, sdn500Times)

    # render the plots

    doPlotPitchRollHeading()
    doPlotVelocity()
    doPlotPosition()
    #doPlotBodyAccel()
    #doPlotBodyRot()
    #doPlotSensorAccel()
    #doPlotSensorRot()

    # show them

    plt.show()

    sys.exit(0)
    
########################################################################
# Read columm headers for the cmigits data
# this is in the first line

def readCmigitsColumnHeaders(filePath):

    colHeaders = []
    colData = {}

    fp = open(filePath, 'r')
    line = fp.readline()
    fp.close()
    
    commentIndex = line.find("#")
    if (commentIndex == 0):
        # header
        colHeaders = line.lstrip("# ").rstrip("\n").split()
        if (options.debug == True):
            print >>sys.stderr, "n colHeaders: ", len(colHeaders)
            print >>sys.stderr, "colHeaders: ", colHeaders
    else:
        print >>sys.stderr, "ERROR - readColumnHeaders"
        print >>sys.stderr, "  First line does not start with #"
        return -1, colHeaders, colData
    
    for index, var in enumerate(colHeaders, start=0):
        colData[var] = []
        
    return 0, colHeaders, colData

########################################################################
# Read columm headers for the sdn500 data
# this is in the first line

def readSdn500ColumnHeaders(filePath):

    colHeaders = []
    colData = {}

    fp = open(filePath, 'r')
    line = fp.readline()
    fp.close()
    
    commentIndex = line.find("#")
    if (commentIndex == 0):
        # header
        colHeaders = line.lstrip("# ").rstrip("\n").split()
        if (options.debug == True):
            print >>sys.stderr, "n colHeaders: ", len(colHeaders)
            print >>sys.stderr, "colHeaders: ", colHeaders
    else:
        print >>sys.stderr, "ERROR - readColumnHeaders"
        print >>sys.stderr, "  First line does not start with #"
        return -1, colHeaders, colData
    
    for index, var in enumerate(colHeaders, start=0):
        colData[var] = []
        
    return 0, colHeaders, colData

########################################################################
# Read in the cmigits data

def readCmigitsInputData(filePath, colHeaders, colData):

    # open file

    fp = open(filePath, 'r')
    lines = fp.readlines()
    lineNum = 0

    # read in a line at a time, set colData
    for line in lines:
        
        commentIndex = line.find("#")
        if (commentIndex >= 0):
            continue

        lineNum = lineNum + 1
        if ((lineNum % int(options.decimation)) != 0):
            continue
            
        # data
        
        data = line.strip().split()
        if (len(data) != len(colHeaders)):
            if (options.debug == True):
                print >>sys.stderr, "skipping line: ", line
                print >>sys.stderr, "  len(data): ", len(data)
            continue;

        for index, var in enumerate(colHeaders, start=0):
            # print >>sys.stderr, "index, data[index]: ", index, ", ", data[index]
            if (var == 'year' or var == 'month' or var == 'day' or \
                var == 'hour' or var == 'minute' or var == 'second' or var == 'microsecond'):
                colData[var].append(int(data[index]))
            else:
                colData[var].append(float(data[index]))

    fp.close()

    # load observation times array

    year = colData['year']
    month = colData['month']
    day = colData['day']
    hour = colData['hour']
    minute = colData['minute']
    sec = colData['second']
    # microsec = colData['microsecond']

    obsTimes = []
    for ii, var in enumerate(year, start=0):
        thisTime = datetime.datetime(year[ii], month[ii], day[ii],
                                     hour[ii], minute[ii], sec[ii],
                                     0)
        obsTimes.append(thisTime)

    return colData, obsTimes

########################################################################
# Read in the sdn500 data

def readSdn500InputData(filePath, colHeaders, colData):

    # open file

    fp = open(filePath, 'r')
    lines = fp.readlines()
    lineNum = 0

    # read in a line at a time, set colData
    for line in lines:
        
        commentIndex = line.find("#")
        if (commentIndex >= 0):
            continue

        lineNum = lineNum + 1
        if ((lineNum % int(options.decimation)) != 0):
            continue
            
        # data
        
        data = line.strip().split()
        if (len(data) != len(colHeaders)):
            if (options.debug == True):
                print >>sys.stderr, "skipping line: ", line
                print >>sys.stderr, "  len(data): ", len(data)
            continue;

        for index, var in enumerate(colHeaders, start=0):
            # print >>sys.stderr, "index, data[index]: ", index, ", ", data[index]
            if (var == 'year' or var == 'month' or var == 'day' or \
                var == 'hour' or var == 'min' or var == 'sec'):
                colData[var].append(int(data[index]))
            else:
                colData[var].append(float(data[index]))

    fp.close()

    # load observation times array

    year = colData['year']
    month = colData['month']
    day = colData['day']
    hour = colData['hour']
    min = colData['min']
    sec = colData['sec']
    # microsec = colData['microsecond']

    obsTimes = []
    for ii, var in enumerate(year, start=0):
        #        print >>sys.stderr, "year, month, day, hour, min, sec: ", \
        #            year[ii], ", ", month[ii], ", ", day[ii], ", ", \
        #            hour[ii], ", ", min[ii], ", ", sec[ii]
        thisTime = datetime.datetime(year[ii], month[ii], day[ii],
                                     hour[ii], min[ii], sec[ii])
        obsTimes.append(thisTime)

    return colData, obsTimes

########################################################################
# Moving average filter

def movingAverage(values, window):

    if (window < 2):
        return values

    weights = np.repeat(1.0, window)/window
    sma = np.convolve(values, weights, 'same')
    return sma

########################################################################
# Set up arrays for plotting

def loadCmigitsDataArrays(cmigitsData, cmigitsTimes):

    filtLen = int(options.filtLen)
    
    # set up arrays

    global cTimes

    cTimes = np.array(cmigitsTimes).astype(datetime.datetime)

    global cLat, cLon, cAlt
    cLat = movingAverage(np.array(cmigitsData["latitude_deg"]).astype(np.double), filtLen)
    cLon = movingAverage(np.array(cmigitsData["longitude_deg"]).astype(np.double), filtLen)
    cAlt = movingAverage(np.array(cmigitsData["altitude_m"]).astype(np.double), filtLen)

    global cVelNorth, cVelEast, cVelUp
    cVelNorth = movingAverage(np.array(cmigitsData["velNorth_m_s-1"]).astype(np.double), filtLen)
    cVelEast = movingAverage(np.array(cmigitsData["velEast_m_s-1"]).astype(np.double), filtLen)
    cVelUp = movingAverage(np.array(cmigitsData["velUp_m_s-1"]).astype(np.double), filtLen)

    global cPitch, cRoll, cHeading
    cPitch = movingAverage(np.array(cmigitsData["pitch_deg"]).astype(np.double), filtLen)
    cRoll = movingAverage(np.array(cmigitsData["roll_deg"]).astype(np.double), filtLen)
    cHeading = movingAverage(np.array(cmigitsData["heading_deg"]).astype(np.double), filtLen)
    
########################################################################
# Set up arrays for plotting

def loadSdn500DataArrays(sdn500Data, sdn500Times):

    filtLen = int(options.filtLen)
    
    # set up arrays

    global sTimes

    sTimes = np.array(sdn500Times).astype(datetime.datetime)

    global sLat, sLon, sAlt
    sLat = movingAverage(np.array(sdn500Data["lat"]).astype(np.double), filtLen)
    sLon = movingAverage(np.array(sdn500Data["lon"]).astype(np.double), filtLen)
    sAlt = movingAverage(np.array(sdn500Data["altKm"]).astype(np.double), filtLen) * 1000.0

    global sVelUp
    sVelUp = movingAverage(np.array(sdn500Data["vertVel"]).astype(np.double), filtLen)

    global sPitch, sRoll, sHeading
    sPitch = movingAverage(np.array(sdn500Data["pitch"]).astype(np.double), filtLen)
    sRoll = movingAverage(np.array(sdn500Data["roll"]).astype(np.double), filtLen)
    sHeading = movingAverage(np.array(sdn500Data["heading"]).astype(np.double), filtLen)
    
########################################################################
# Plot pitch, roll and heading data

def doPlotPitchRollHeading():

    # set up plots

    widthIn = float(options.mainWidthMm) / 25.4
    htIn = float(options.mainHeightMm) / 25.4
    
    global figNum
    fig = plt.figure(figNum, (widthIn, htIn))
    figNum = figNum + 1
    
    ax1 = fig.add_subplot(3,1,1,xmargin=0.0)
    ax2 = fig.add_subplot(3,1,2,xmargin=0.0)
    ax3 = fig.add_subplot(3,1,3,xmargin=0.0)
    
    ax1.plot(cTimes, cPitch, label='Cmigits Pitch(deg)', color='red', linewidth=1)
    ax1.plot(sTimes, sPitch, label='Sdn500  Pitch(deg)', color='blue', linewidth=1)

    ax2.plot(cTimes, cRoll, label='Cmigits Roll(deg)', color='red', linewidth=1)
    ax2.plot(sTimes, sRoll, label='Sdn500  Roll(deg)', color='blue', linewidth=1)

    ax3.plot(cTimes, cHeading, label='Cmigits Heading(deg)', color='red', linewidth=1)
    ax3.plot(sTimes, sHeading, label='Sdn500  Heading(deg)', color='blue', linewidth=1)

    configTimeAxis(ax1, -9999, -9999, "Pitch(deg)", 'upper right')
    configTimeAxis(ax2, -9999, -9999, "Roll(deg)", 'upper right')
    configTimeAxis(ax3, -9999, -9999, "Heading(deg)", 'upper right')

    fig.autofmt_xdate()
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.08, left=0.06, right=0.97, top=0.96)
    
    # title name

    fig.suptitle("PITCH/ROLL/HDG - " + str(startTime) + " to " + str(endTime))

    return

########################################################################
# Plot velocities

def doPlotVelocity():

    # set up plots

    widthIn = float(options.mainWidthMm) / 25.4
    htIn = float(options.mainHeightMm) / 25.4
    
    global figNum
    fig = plt.figure(figNum, (widthIn, htIn))
    figNum = figNum + 1
    
    ax1 = fig.add_subplot(3,1,1,xmargin=0.0)
    ax2 = fig.add_subplot(3,1,2,xmargin=0.0)
    ax3 = fig.add_subplot(3,1,3,xmargin=0.0)
    
    ax1.plot(cTimes, cVelNorth, label='Cmigits VelNorth(m/s)', color='red', linewidth=1)

    ax2.plot(cTimes, cVelEast, label='Cmigits VelEast(m/s)', color='red', linewidth=1)

    ax3.plot(cTimes, cVelUp, label='Cmigits VelUp(m/s)', color='red', linewidth=1)
    ax3.plot(sTimes, sVelUp, label='Sdn500 VelUp(m/s)', color='blue', linewidth=1)

    configTimeAxis(ax1, -9999, -9999, "VelNorth", 'upper right')
    configTimeAxis(ax2, -9999, -9999, "VelEast", 'upper right')
    configTimeAxis(ax3, -9999, -9999, "VelUp", 'upper right')

    fig.autofmt_xdate()
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.08, left=0.06, right=0.97, top=0.96)
    
    # title name

    fig.suptitle("VELOCITIES - " + str(startTime) + " to " + str(endTime))

    return

########################################################################
# Plot position

def doPlotPosition():

    # set up plots

    widthIn = float(options.mainWidthMm) / 25.4
    htIn = float(options.mainHeightMm) / 25.4
    
    global figNum
    fig = plt.figure(figNum, (widthIn, htIn))
    figNum = figNum + 1
    
    ax1 = fig.add_subplot(3,1,1,xmargin=0.0)
    ax2 = fig.add_subplot(3,1,2,xmargin=0.0)
    ax3 = fig.add_subplot(3,1,3,xmargin=0.0)
    
    ax1.plot(cTimes, cLat, label='cmigits latitude', color='red', linewidth=1)
    ax1.plot(sTimes, sLat, label='sdn500 latitude', color='blue', linewidth=1)

    ax2.plot(cTimes, cLon, label='cmigits longitude', color='red', linewidth=1)
    ax2.plot(sTimes, sLon, label='sdn500 longitude', color='blue', linewidth=1)

    ax3.plot(cTimes, cAlt, label='cmigits altitude(m)', color='red', linewidth=1)
    ax3.plot(sTimes, sAlt, label='sdn500 altitude(m)', color='blue', linewidth=1)
    
    configTimeAxis(ax1, -9999, -9999, "latitude", 'upper right')
    configTimeAxis(ax2, -9999, -9999, "longitude", 'upper right')
    #configTimeAxis(ax1, 39, 43, "latitude", 'upper right')
    #configTimeAxis(ax2, -106, -100, "longitude", 'upper right')
    configTimeAxis(ax3, -9999, -9999, "altitude", 'upper right')

    fig.autofmt_xdate()
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.08, left=0.06, right=0.97, top=0.96)
    
    # title name

    fig.suptitle("FOG LAT/LON/ALT POSITION - " + str(startTime) + " to " + str(endTime))

    return

########################################################################
# Configure axes, legends etc

def configTimeAxis(ax, miny, maxy, ylabel, legendLoc):
    
    legend = ax.legend(loc=legendLoc, ncol=8)
    for label in legend.get_texts():
        label.set_fontsize('x-small')
        ax.set_xlabel("Time")
    ax.set_ylabel(ylabel)
    ax.grid(True)
    if (miny > -9990 and maxy > -9990):
        ax.set_ylim([miny, maxy])
    hfmt = dates.DateFormatter('%H:%M:%S')
    ax.xaxis.set_major_locator(dates.AutoDateLocator())
    ax.xaxis.set_major_formatter(hfmt)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(8) 

    ax.set_xlim(startTime, endTime)

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

