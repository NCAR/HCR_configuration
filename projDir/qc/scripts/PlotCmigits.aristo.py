#!/usr/bin/env python

#===========================================================================
#
# Produce plots for SPATIAL FOG data for ARISTO
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

    appName = 'PlotSpatialFog.py'
    projDir = os.environ['PROJ_DIR']
    dataDir = os.environ['DATA_DIR']
    global dataFilePath
    dataFilePath = os.path.join(dataDir, 'fog/SpatialFog_RF04.txt')

    # parse the command line

    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option('--debug',
                      dest='debug', default=False,
                      action="store_true",
                      help='Set debugging on')
    parser.add_option('--file',
                      dest='dataFilePath',
                      default=dataFilePath,
                      help='File path for comparison results')
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
                      default=21,
                      help='Len of moving mean filter')
    parser.add_option('--decimation',
                      dest='decimation',
                      default=20,
                      help='Decimation factor when reading in data')
    parser.add_option('--start',
                      dest='startTime',
                      default='2016 10 10 16 00 00',
                      help='Start time for XY plot')
    parser.add_option('--end',
                      dest='endTime',
                      default='2016 10 10 17 45 00',
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
        print >>sys.stderr, "  dataFilePath: ", options.dataFilePath
        print >>sys.stderr, "  startTime: ", startTime
        print >>sys.stderr, "  endTime: ", endTime

    # read in column headers for bias results

    iret, compHdrs, compData = readColumnHeaders(options.dataFilePath)
    if (iret != 0):
        sys.exit(-1)

    # read in data for comp results

    compData, compTimes = readInputData(options.dataFilePath, compHdrs, compData)

    # load up the data arrays

    loadDataArrays(compData, compTimes)

    # render the plots

    doPlotPitchRollHeading()
    #doPlotVelocity()
    doPlotPosition()
    #doPlotBodyAccel()
    #doPlotBodyRot()
    #doPlotSensorAccel()
    #doPlotSensorRot()

    # show them

    plt.show()

    sys.exit(0)
    
########################################################################
# Read columm headers for the data
# this is in the first line

def readColumnHeaders(filePath):

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
# Read in the data

def readInputData(filePath, colHeaders, colData):

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
# Moving average filter

def movingAverage(values, window):

    if (window < 2):
        return values

    weights = np.repeat(1.0, window)/window
    sma = np.convolve(values, weights, 'same')
    return sma

########################################################################
# Set up arrays for plotting

def loadDataArrays(compData, compTimes):

    filtLen = int(options.filtLen)
    
    # set up arrays

    global ctimes

    ctimes = np.array(compTimes).astype(datetime.datetime)

    global lat, lon, alt
    lat = movingAverage(np.array(compData["latitude_deg"]).astype(np.double), filtLen)
    lon = movingAverage(np.array(compData["longitude_deg"]).astype(np.double), filtLen)
    alt = movingAverage(np.array(compData["altitude_m"]).astype(np.double), filtLen)

    global velNorth, velEast, velUp
    velNorth = movingAverage(np.array(compData["velNorth_m_s-1"]).astype(np.double), filtLen)
    velEast = movingAverage(np.array(compData["velEast_m_s-1"]).astype(np.double), filtLen)
    velUp = movingAverage(np.array(compData["velUp_m_s-1"]).astype(np.double), filtLen)

    global pitch, roll, heading
    pitch = movingAverage(np.array(compData["pitch_deg"]).astype(np.double), filtLen)
    roll = movingAverage(np.array(compData["roll_deg"]).astype(np.double), filtLen)
    heading = movingAverage(np.array(compData["heading_deg"]).astype(np.double), filtLen)
    
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
    
    ax1.plot(ctimes, pitch, label='Pitch(deg)', color='green', linewidth=1)
    ax2.plot(ctimes, roll, label='Roll(deg)', color='red', linewidth=1)
    ax3.plot(ctimes, heading, label='Heading(deg)', color='blue', linewidth=1)

    configTimeAxis(ax1, -9999, -9999, "Pitch(deg)", 'upper right')
    configTimeAxis(ax2, -9999, -9999, "Roll(deg)", 'upper right')
    configTimeAxis(ax3, -9999, -9999, "Heading(deg)", 'upper right')

    fig.autofmt_xdate()
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.08, left=0.06, right=0.97, top=0.96)
    
    # title name

    fig.suptitle("FOG PITCH/ROLL/HDG - " + str(startTime) + " to " + str(endTime))

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
    
    ax1.plot(ctimes, velNorth, label='VelNorth(m/s)', color='green', linewidth=1)
    ax2.plot(ctimes, velEast, label='VelEast(m/s)', color='red', linewidth=1)
    ax3.plot(ctimes, velUp, label='VelUp(m/s)', color='blue', linewidth=1)

    configTimeAxis(ax1, -9999, -9999, "VelNorth", 'upper right')
    configTimeAxis(ax2, -9999, -9999, "VelEast", 'upper right')
    configTimeAxis(ax3, -9999, -9999, "VelUp", 'upper right')

    fig.autofmt_xdate()
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.08, left=0.06, right=0.97, top=0.96)
    
    # title name

    fig.suptitle("FOG VELOCITIES - " + str(startTime) + " to " + str(endTime))

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
    
    ax1.plot(ctimes, lat, label='latitude', color='green', linewidth=1)
    ax2.plot(ctimes, lon, label='longitude', color='red', linewidth=1)
    ax3.plot(ctimes, alt, label='altitude(m)', color='blue', linewidth=1)

    #configTimeAxis(ax1, -9999, -9999, "latitude", 'upper right')
    #configTimeAxis(ax2, -9999, -9999, "longitude", 'upper right')
    configTimeAxis(ax1, 39, 43, "latitude", 'upper right')
    configTimeAxis(ax2, -106, -100, "longitude", 'upper right')
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

