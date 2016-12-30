#!/usr/bin/env python

#===========================================================================
#
# Produce plots for GEOREF differences
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

    appName = 'PlotGeorefDiffs.fog.py'
    global compFilePath
    compFilePath = '/tmp/' + appName + '.' + str(os.getpid()) + '.txt'

    # parse the command line

    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option('--debug',
                      dest='debug', default=False,
                      action="store_true",
                      help='Set debugging on')
    parser.add_option('--file',
                      dest='compFilePath',
                      default=compFilePath,
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
    parser.add_option('--start',
                      dest='startTime',
                      default='2016 09 19 16 50 00',
                      help='Start time for XY plot')
    parser.add_option('--widthSecs',
                      dest='widthSecs',
                      default=6000,
                      help='Width of plot in seconds')
    
    (options, args) = parser.parse_args()
    
    # time limits

    year, month, day, hour, minute, sec = options.startTime.split()
    startTime = datetime.datetime(int(year), int(month), int(day),
                                  int(hour), int(minute), int(sec))
    widthTimeDelta = datetime.timedelta(0, int(options.widthSecs))
    endTime = startTime + widthTimeDelta

    if (options.debug == True):
        print >>sys.stderr, "  compFilePath: ", options.compFilePath
        print >>sys.stderr, "  startTime: ", startTime
        print >>sys.stderr, "  endTime: ", endTime

    # read in column headers for bias results

    iret, compHdrs, compData = readColumnHeaders(options.compFilePath)
    if (iret != 0):
        sys.exit(-1)

    # read in data for comp results

    compData, compTimes = readInputData(options.compFilePath, compHdrs, compData)

    # load up the data arrays

    loadDataArrays(compData, compTimes)

    # render the plots

    doPlotPitch()
    doPlotRoll()
    doPlotHeading()
    #doPlotTrack()
    doPlotDrift()
    doPlotVertVel()

    #doPlotAltKm()
    #doPlotLat()
    #doPlotLon()

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

    # read in a line at a time, set colData
    for line in lines:
        
        commentIndex = line.find("#")
        if (commentIndex >= 0):
            continue
            
        # data
        
        data = line.strip().split()
        if (len(data) != len(colHeaders)):
            if (options.debug == True):
                print >>sys.stderr, "skipping line: ", line
            continue;

        for index, var in enumerate(colHeaders, start=0):
            # print >>sys.stderr, "index, data[index]: ", index, ", ", data[index]
            if (var == 'count' or var == 'year' or var == 'month' or var == 'day' or \
                var == 'hour' or var == 'min' or var == 'sec' or \
                var == 'unix_time'):
                colData[var].append(int(data[index]))
            else:
                colData[var].append(float(data[index]))

    fp.close()

    # load observation times array

    year = colData['year']
    month = colData['month']
    day = colData['day']
    hour = colData['hour']
    minute = colData['min']
    sec = colData['sec']

    obsTimes = []
    for ii, var in enumerate(year, start=0):
        thisTime = datetime.datetime(year[ii], month[ii], day[ii],
                                     hour[ii], minute[ii], sec[ii])
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

    global pitch, pitchSec, pitchDiff
    pitch = movingAverage(np.array(compData["pitch"]).astype(np.double), filtLen)
    pitchSec = movingAverage(np.array(compData["pitchSec"]).astype(np.double), filtLen)
    pitchDiff = movingAverage(np.array(compData["pitchDiff"]).astype(np.double), filtLen)

    global roll, rollSec, rollDiff
    roll = movingAverage(np.array(compData["roll"]).astype(np.double), filtLen)
    rollSec = movingAverage(np.array(compData["rollSec"]).astype(np.double), filtLen)
    rollDiff = movingAverage(np.array(compData["rollDiff"]).astype(np.double), filtLen)
    
    global heading, headingSec, headingDiff
    heading = movingAverage(np.array(compData["heading"]).astype(np.double), filtLen)
    headingSec = movingAverage(np.array(compData["headingSec"]).astype(np.double), filtLen)
    headingDiff = movingAverage(np.array(compData["hdgDiff"]).astype(np.double), filtLen)
    
    global drift, driftSec, driftDiff
    drift = movingAverage(np.array(compData["drift"]).astype(np.double), filtLen)
    driftSec = movingAverage(np.array(compData["driftSec"]).astype(np.double), filtLen)
    driftDiff = movingAverage(np.array(compData["driftDiff"]).astype(np.double), filtLen)

    global track, trackSec, trackDiff
    track = np.array(compData["track"]).astype(np.double)
    trackSec = np.array(compData["trackSec"]).astype(np.double)
    trackDiff = np.array(compData["trackDiff"]).astype(np.double)

    global altKm, altKmSec, altKmDiff
    altKm = np.array(compData["altKm"]).astype(np.double)
    altKmSec = np.array(compData["altKmSec"]).astype(np.double)
    altKmDiff = np.array(compData["altKmDiff"]).astype(np.double)

    global vertVel, vertVelSec, vertVelDiff
    vertVel = np.array(compData["vertVel"]).astype(np.double)
    vertVelSec = np.array(compData["vertVelSec"]).astype(np.double)
    vertVelDiff = np.array(compData["vertVelDiff"]).astype(np.double)

    global lat, latSec, latDiff
    lat = np.array(compData["lat"]).astype(np.double)
    latSec = np.array(compData["latSec"]).astype(np.double)
    latDiff = np.array(compData["latDiff"]).astype(np.double)

    global lon, lonSec, lonDiff
    lon = np.array(compData["lon"]).astype(np.double)
    lonSec = np.array(compData["lonSec"]).astype(np.double)
    lonDiff = np.array(compData["lonDiff"]).astype(np.double)

########################################################################
# Plot pitch data

def doPlotPitch():

    # set up plots

    widthIn = float(options.mainWidthMm) / 25.4
    htIn = float(options.mainHeightMm) / 25.4
    
    global figNum
    fig = plt.figure(figNum, (widthIn, htIn))
    figNum = figNum + 1
    
    ax1 = fig.add_subplot(2,1,1,xmargin=0.0)
    ax2 = fig.add_subplot(2,1,2,xmargin=0.0)
    
    ax1.plot(ctimes, pitch, label='Pitch-INS1', color='green', linewidth=1)
    ax1.plot(ctimes, pitchSec, label='Pitch-FOG', color='red', linewidth=1)
    
    ax2.plot(ctimes, pitchDiff, \
             label='pitchDiff', color='blue', linewidth=1)
    
    configTimeAxis(ax1, -9999, -9999, "Pitch", 'upper right')
    configTimeAxis(ax2, -9999, -9999, "PitchDiff", 'upper right')
    #configTimeAxis(ax1, -2, 5, "Pitch", 'upper right')
    #configTimeAxis(ax2, -2, 4, "PitchDiff", 'upper right')

    fig.autofmt_xdate()
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.08, left=0.06, right=0.97, top=0.96)
    
    # title name

    fig.suptitle("PITCH - " + str(startTime) + " to " + str(endTime))

    return

########################################################################
# Plot roll data

def doPlotRoll():

    # set up plots

    widthIn = float(options.mainWidthMm) / 25.4
    htIn = float(options.mainHeightMm) / 25.4
    
    global figNum
    fig = plt.figure(figNum, (widthIn, htIn))
    figNum = figNum + 1
    
    ax1 = fig.add_subplot(2,1,1,xmargin=0.0)
    ax2 = fig.add_subplot(2,1,2,xmargin=0.0)
    
    ax1.plot(ctimes, roll, label='Roll-INS1', color='green', linewidth=1)
    ax1.plot(ctimes, rollSec, label='Roll-FOG', color='red', linewidth=1)
    
    ax2.plot(ctimes, rollDiff, \
             label='rollDiff', color='blue', linewidth=1)
    
    configTimeAxis(ax1, -9999, -9999, "Roll", 'upper right')
    configTimeAxis(ax2, -9999, -9999, "RollDiff", 'upper right')

    fig.autofmt_xdate()
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.08, left=0.06, right=0.97, top=0.96)
    
    # title name

    fig.suptitle("ROLL - " + str(startTime) + " to " + str(endTime))

    return

########################################################################
# Plot heading data

def doPlotHeading():

    # set up plots

    widthIn = float(options.mainWidthMm) / 25.4
    htIn = float(options.mainHeightMm) / 25.4
    
    global figNum
    fig = plt.figure(figNum, (widthIn, htIn))
    figNum = figNum + 1
    
    ax1 = fig.add_subplot(2,1,1,xmargin=0.0)
    ax2 = fig.add_subplot(2,1,2,xmargin=0.0)
    
    ax1.plot(ctimes, heading, label='Heading-INS1', color='green', linewidth=1)
    ax1.plot(ctimes, headingSec, label='Heading-FOG', color='red', linewidth=1)
    
    ax2.plot(ctimes, headingDiff, \
             label='headingDiff', color='blue', linewidth=1)
    
    configTimeAxis(ax1, -9999, -9999, "Heading", 'upper right')
    #configTimeAxis(ax2, -9999, -9999, "HeadingDiff", 'upper right')
    configTimeAxis(ax2, -5, 5, "HeadingDiff", 'upper right')

    fig.autofmt_xdate()
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.08, left=0.06, right=0.97, top=0.96)
    
    # title name

    fig.suptitle("HEADING - " + str(startTime) + " to " + str(endTime))

    return

########################################################################
# Plot track data

def doPlotTrack():

    # set up plots

    widthIn = float(options.mainWidthMm) / 25.4
    htIn = float(options.mainHeightMm) / 25.4
    
    global figNum
    fig = plt.figure(figNum, (widthIn, htIn))
    figNum = figNum + 1
    
    ax1 = fig.add_subplot(2,1,1,xmargin=0.0)
    ax2 = fig.add_subplot(2,1,2,xmargin=0.0)
    
    ax1.plot(ctimes, track, label='Track-INS1', color='green', linewidth=1)
    ax1.plot(ctimes, trackSec, label='Track-FOG', color='red', linewidth=1)
    
    ax2.plot(ctimes, trackDiff, \
             label='trackDiff', color='blue', linewidth=1)
    
    configTimeAxis(ax1, -9999, -9999, "Track", 'upper right')
    configTimeAxis(ax2, -9999, -9999, "TrackDiff", 'upper right')

    fig.autofmt_xdate()
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.08, left=0.06, right=0.97, top=0.96)
    
    # title name

    fig.suptitle("TRACK - " + str(startTime) + " to " + str(endTime))

    return

########################################################################
# Plot drift data

def doPlotDrift():

    # set up plots

    widthIn = float(options.mainWidthMm) / 25.4
    htIn = float(options.mainHeightMm) / 25.4
    
    global figNum
    fig = plt.figure(figNum, (widthIn, htIn))
    figNum = figNum + 1
    
    ax1 = fig.add_subplot(2,1,1,xmargin=0.0)
    ax2 = fig.add_subplot(2,1,2,xmargin=0.0)
    
    ax1.plot(ctimes, drift, label='Drift-INS1', color='green', linewidth=1)
    ax1.plot(ctimes, driftSec, label='Drift-FOG', color='red', linewidth=1)
    
    ax2.plot(ctimes, driftDiff, \
             label='driftDiff', color='blue', linewidth=1)
    
    configTimeAxis(ax1, -9999, -9999, "Drift", 'upper right')
    configTimeAxis(ax2, -9999, -9999, "DriftDiff", 'upper right')

    fig.autofmt_xdate()
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.08, left=0.06, right=0.97, top=0.96)
    
    # title name

    fig.suptitle("DRIFT - " + str(startTime) + " to " + str(endTime))

    return

########################################################################
# Plot vertVel data

def doPlotVertVel():

    # set up plots

    widthIn = float(options.mainWidthMm) / 25.4
    htIn = float(options.mainHeightMm) / 25.4
    
    global figNum
    fig = plt.figure(figNum, (widthIn, htIn))
    figNum = figNum + 1
    
    ax1 = fig.add_subplot(2,1,1,xmargin=0.0)
    ax2 = fig.add_subplot(2,1,2,xmargin=0.0)
    
    ax1.plot(ctimes, vertVel, label='VertVel-INS1', color='green', linewidth=1)
    ax1.plot(ctimes, vertVelSec, label='VertVel-FOG', color='red', linewidth=1)
    
    ax2.plot(ctimes, vertVelDiff, \
             label='vertVelDiff', color='blue', linewidth=1)
    
    configTimeAxis(ax1, -9999, -9999, "VertVel", 'upper right')
    configTimeAxis(ax2, -9999, -9999, "VertVelDiff", 'upper right')

    fig.autofmt_xdate()
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.08, left=0.06, right=0.97, top=0.96)
    
    # title name

    fig.suptitle("VERTVEL - " + str(startTime) + " to " + str(endTime))

    return

########################################################################
# Plot altKm data

def doPlotAltKm():

    # set up plots

    widthIn = float(options.mainWidthMm) / 25.4
    htIn = float(options.mainHeightMm) / 25.4
    
    global figNum
    fig = plt.figure(figNum, (widthIn, htIn))
    figNum = figNum + 1
    
    ax1 = fig.add_subplot(2,1,1,xmargin=0.0)
    ax2 = fig.add_subplot(2,1,2,xmargin=0.0)
    
    ax1.plot(ctimes, altKm, label='AltKm-INS1', color='green', linewidth=1)
    ax1.plot(ctimes, altKmSec, label='AltKm-FOG', color='red', linewidth=1)
    
    ax2.plot(ctimes, altKmDiff, \
             label='altKmDiff', color='blue', linewidth=1)
    
    configTimeAxis(ax1, -9999, -9999, "AltKm", 'upper right')
    configTimeAxis(ax2, -9999, -9999, "AltKmDiff", 'upper right')

    fig.autofmt_xdate()
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.08, left=0.06, right=0.97, top=0.96)
    
    # title name

    fig.suptitle("ALTKM - " + str(startTime) + " to " + str(endTime))

    return

########################################################################
# Plot lat data
 
def doPlotLat():

    # set up plots

    widthIn = float(options.mainWidthMm) / 25.4
    htIn = float(options.mainHeightMm) / 25.4
    
    global figNum
    fig = plt.figure(figNum, (widthIn, htIn))
    figNum = figNum + 1
    
    ax1 = fig.add_subplot(2,1,1,xmargin=0.0)
    ax2 = fig.add_subplot(2,1,2,xmargin=0.0)
    
    ax1.plot(ctimes, lat, label='Lat-INS1', color='green', linewidth=1)
    ax1.plot(ctimes, latSec, label='Lat-FOG', color='red', linewidth=1)
    
    ax2.plot(ctimes, latDiff, \
             label='latDiff', color='blue', linewidth=1)
    
    configTimeAxis(ax1, -9999, -9999, "Lat", 'upper right')
    configTimeAxis(ax2, -9999, -9999, "LatDiff", 'upper right')

    fig.autofmt_xdate()
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.08, left=0.06, right=0.97, top=0.96)
    
    # title name

    fig.suptitle("LATITUDE - " + str(startTime) + " to " + str(endTime))

    return

########################################################################
# Plot lon data
 
def doPlotLon():

    # set up plots

    widthIn = float(options.mainWidthMm) / 25.4
    htIn = float(options.mainHeightMm) / 25.4
    
    global figNum
    fig = plt.figure(figNum, (widthIn, htIn))
    figNum = figNum + 1
    
    ax1 = fig.add_subplot(2,1,1,xmargin=0.0)
    ax2 = fig.add_subplot(2,1,2,xmargin=0.0)
    
    ax1.plot(ctimes, lon, label='Lon-INS1', color='green', linewidth=1)
    ax1.plot(ctimes, lonSec, label='Lon-FOG', color='red', linewidth=1)
    
    ax2.plot(ctimes, lonDiff, \
             label='lonDiff', color='blue', linewidth=1)
    
    configTimeAxis(ax1, -9999, -9999, "Lon", 'upper right')
    configTimeAxis(ax2, -9999, -9999, "LonDiff", 'upper right')

    fig.autofmt_xdate()
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.08, left=0.06, right=0.97, top=0.96)
    
    # title name

    fig.suptitle("LONGITUDE - " + str(startTime) + " to " + str(endTime))

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

