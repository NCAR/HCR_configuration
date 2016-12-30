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

    appName = 'PlotGeorefDiffs.iwg1_vs_10hz.py'
    global compFilePath
    compFilePath = '/tmp/' + appName + '.' + str(os.getpid()) + '.txt'

    # parse the command line

    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option('--debug',
                      dest='debug', default=False,
                      action="store_true",
                      help='Set debugging on')
    parser.add_option('--verbose',
                      dest='verbose', default=False,
                      action="store_true",
                      help='Set verbose debugging on')
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
                      default=11,
                      help='Len of moving mean filter')
    parser.add_option('--start',
                      dest='startTime',
                      default='2015 07 27 20 00 00',
                      help='Start time for XY plot')
    parser.add_option('--widthSecs',
                      dest='widthSecs',
                      default=1800,
                      help='Width of plot in seconds')
    parser.add_option('--realtime',
                      dest='realtime',
                      default=False,
                      action="store_true",
                      help='Run in realtime mode. End time is now')
    
    (options, args) = parser.parse_args()
    
    if (options.verbose == True):
        options.debug = True

    # time limits

    if (options.realtime):
        now = datetime.datetime.utcnow()
        widthTimeDelta = datetime.timedelta(0, int(options.widthSecs))
        start = now - widthTimeDelta
        startTime = datetime.datetime(start.year, start.month, start.day, \
                                      start.hour, start.minute, start.second)
        endTime = datetime.datetime(now.year, now.month, now.day, \
                                    now.hour, now.minute, now.second)
    else:
        year, month, day, hour, minute, sec = options.startTime.split()
        startTime = datetime.datetime(int(year), int(month), int(day),
                                      int(hour), int(minute), int(sec))
        widthTimeDelta = datetime.timedelta(0, int(options.widthSecs))
        endTime = startTime + widthTimeDelta

    if (options.debug == True):
        print >>sys.stderr, "Running %prog"
        print >>sys.stderr, "  compFilePath: ", options.compFilePath
        if (options.realtime):
            print >>sys.stderr, "  end time is now: ", now
        print >>sys.stderr, "  startTime: ", startTime
        print >>sys.stderr, "  endTime: ", endTime

    # create the comparison text file

    createCompFile()

    # read in column headers for bias results

    iret, compHdrs, compData = readColumnHeaders(options.compFilePath)
    if (iret != 0):
        sys.exit(-1)

    # read in data for comp results

    compData, compTimes = readInputData(options.compFilePath, compHdrs, compData)

    # load up the data arrays

    loadDataArrays(compData, compTimes)

    # render the plots
    
    #doPlotOverview()
    doPlotPitchRoll()
    doPlotDiffs()
    #doPlotEstPitchDiff()
    #doPlotRadarAngles()

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

    global temp, tempHcrFog, tempTailcone
    global pressure, rh, density, altPres, altGps
    global vertVel, gvMassKg, gvMass10000Kg
    global aoa, aoaSm, aoa2, ias, tas

    temp = np.array(compData["temp"]).astype(np.double)
    tempHcrFog = np.array(compData["tempSec"]).astype(np.double)
    tempTailcone = np.array(compData["custom0Sec"]).astype(np.double)
    tempNose = np.array(compData["custom6"]).astype(np.double)
    tempAtf1 = np.array(compData["custom7"]).astype(np.double)
    pressure = np.array(compData["pressure"]).astype(np.double)
    rh = np.array(compData["rh"]).astype(np.double)
    density = np.array(compData["density"]).astype(np.double)
    altPres = np.array(compData["altPresM"]).astype(np.double)
    altGps = np.array(compData["altGpsM"]).astype(np.double)
    vertVel = np.array(compData["vertVel"]).astype(np.double)
    gvMassKg = np.array(compData["weightKg"]).astype(np.double)
    gvMass10000Kg = np.array(compData["weightKg"]).astype(np.double) / 10000.0
    aoa = np.array(compData["aoa"]).astype(np.double)
    aoaSm = movingAverage(aoa, filtLen)
    aoa2 = aoaSm / 2.0
    ias = np.array(compData["ias"]).astype(np.double)
    tas = np.array(compData["tas"]).astype(np.double)

    global temp30Offset, temp30OffsetSm, temp30OffsetBy10
    temp30Offset = 30.0 - tempHcrFog
    temp30OffsetSm = movingAverage(temp30Offset, filtLen)
    temp30OffsetBy10 = temp30OffsetSm / 10.0

    global accelNorm, accelNormSm, accelLat, accelLatSm, accelLon, accelLonSm

    accelNorm = np.array(compData["accelNorm"]).astype(np.double)
    accelNormSm = movingAverage(accelNorm, filtLen) + 1.0
    accelLat = np.array(compData["accelLat"]).astype(np.double)
    accelLatSm = movingAverage(accelLat, filtLen)
    accelLon = np.array(compData["accelLon"]).astype(np.double)
    accelLonSm = movingAverage(accelLon, filtLen)

    global altDiff, pitchDiff, rollDiff, trackDiff, hdgDiff, vVelDiff
    global pitchDiffSm, rollDiffSm, trackDiffSm, hdgDiffSm, vVelDiffSm
    global driftDiff, driftDiffSm

    altDiff = np.array(compData["altKmDiff"]).astype(np.double)
    altDiffSm = movingAverage(altDiff, filtLen)

    pitchDiff = np.array(compData["pitchDiff"]).astype(np.double)
    pitchDiffSm = movingAverage(pitchDiff, filtLen)

    rollDiff = np.array(compData["rollDiff"]).astype(np.double)
    rollDiffSm = movingAverage(rollDiff, filtLen)

    trackDiff = np.array(compData["trackDiff"]).astype(np.double)
    trackDiffSm = movingAverage(trackDiff, filtLen)

    hdgDiff = np.array(compData["hdgDiff"]).astype(np.double)
    hdgDiffSm = movingAverage(hdgDiff, filtLen)

    driftDiff = np.array(compData["driftDiff"]).astype(np.double)
    driftDiffSm = movingAverage(driftDiff, filtLen)

    vVelDiff = np.array(compData["vertVelDiff"]).astype(np.double)
    vVelDiffSm = movingAverage(vVelDiff, filtLen)
    
    global pitch, pitch2, pitch3, pitchSm, pitch2Sm, pitch3Sm
    global pitchDiff2, pitchDiff2Sm, pitchDiff3, pitchDiff3Sm
    pitch = np.array(compData["pitch"]).astype(np.double)
    pitchSm = movingAverage(pitch, filtLen)
    pitch2 = np.array(compData["custom0"]).astype(np.double)
    pitch2Sm = movingAverage(pitch2, filtLen)
    pitch3 = np.array(compData["custom1"]).astype(np.double)
    pitch3Sm = movingAverage(pitch3, filtLen)
    pitchDiff2 = pitch - pitch2
    pitchDiff2Sm = movingAverage(pitchDiff2, filtLen)
    pitchDiff3 = pitch - pitch3
    pitchDiff3Sm = movingAverage(pitchDiff3, filtLen)

    global pitchCm, pitchCmSm, pitchDiffCm2, pitchDiffCm2Sm
    pitchCm = np.array(compData["pitchSec"]).astype(np.double)
    pitchCmSm = movingAverage(pitchCm, filtLen)
    pitchDiffCm2 = pitch2 - pitchCm - 0.15
    pitchDiffCm2Sm = movingAverage(pitchDiffCm2, filtLen)

    global roll, roll2, roll3, rollSm, roll2Sm, roll3Sm
    global rollDiff2, rollDiff2Sm, rollDiff3, rollDiff3Sm
    roll = np.array(compData["roll"]).astype(np.double)
    rollSm = movingAverage(roll, filtLen)
    roll2 = np.array(compData["custom2"]).astype(np.double)
    roll2Sm = movingAverage(roll2, filtLen)
    roll3 = np.array(compData["custom3"]).astype(np.double)
    roll3Sm = movingAverage(roll3, filtLen)
    rollDiff2 = roll - roll2
    rollDiff2Sm = movingAverage(rollDiff2, filtLen)
    rollDiff3 = roll - roll3
    rollDiff3Sm = movingAverage(rollDiff3, filtLen)

    global rollCm, rollCmSm
    rollCm = np.array(compData["rollSec"]).astype(np.double)
    rollCmSm = movingAverage(rollCm, filtLen)

    global drift, drift2, drift3, driftDiff2, driftDiff2Sm, driftDiff3, driftDiff3Sm
    drift = np.array(compData["drift"]).astype(np.double)
    drift2 = np.array(compData["custom4"]).astype(np.double)
    drift3 = np.array(compData["custom5"]).astype(np.double)
    driftDiff2 = drift - drift2
    driftDiff2Sm = movingAverage(driftDiff2, filtLen)
    driftDiff3 = drift - drift3
    driftDiff3Sm = movingAverage(driftDiff3, filtLen)

    global surfaceVel, surfaceVelSm
    surfaceVel = np.array(compData["custom1Sec"]).astype(np.double)
    surfaceVelSm = movingAverage(surfaceVel, filtLen * 5)

    global azimuth, elevation, rotation, tilt
    azimuth = np.array(compData["custom2Sec"]).astype(np.double)
    elevation = np.array(compData["custom3Sec"]).astype(np.double)
    rotation = np.array(compData["custom4Sec"]).astype(np.double)
    tilt = np.array(compData["custom5Sec"]).astype(np.double)

    global elevErr, elevErrSm, tiltErr, tiltErrSm
    elevErr = elevation
    tiltErr = tilt
    for index, elev in enumerate(elevation):
        if (elev < -85.0):
            elevErr[index] = -90.0 - elevation[index]
            tiltErr[index] = pitchCm[index] + tilt[index]
        if (elev > 85.0):
            elevErr[index] = 90.0 - elevation[index]
            tiltErr[index] = pitchCm[index] - tilt[index]
    elevErrSm = movingAverage(elevErr, filtLen)
    tiltErrSm = movingAverage(tiltErr, filtLen)

    global estPitchDiff, estPitchDiffSm
    estPitchDiff = pitchDiff
    for index, vel in enumerate(surfaceVel):
        tasVal = tas[index]
        if (np.isfinite(tasVal) and np.isfinite(vel) and math.fabs(vel) < 1000):
            angErr = math.degrees(math.asin(vel / tasVal))
            estPitchDiff[index] = pitchDiff[index] + angErr
        else:
            estPitchDiff[index] = float('nan')

    estPitchDiffSm = movingAverage(estPitchDiff, filtLen)

########################################################################
# Plot aircraft variables

def doPlotOverview():

    # set up plots

    widthIn = float(options.mainWidthMm) / 25.4
    htIn = float(options.mainHeightMm) / 25.4

    global figNum
    fig = plt.figure(figNum, (widthIn, htIn))
    figNum = figNum + 1
    
    ax1 = fig.add_subplot(4,1,1,xmargin=0.0)
    ax2 = fig.add_subplot(4,1,2,xmargin=0.0)
    ax3 = fig.add_subplot(4,1,3,xmargin=0.0)
    ax4 = fig.add_subplot(4,1,4,xmargin=0.0)
    
    ax1.plot(ctimes, altPres, \
             label='Pressure altitude (m)', color='blue', linewidth=1)
    
    ax2.plot(ctimes, temp, label='Air temp (C)', color='red', linewidth=1)
    ax2.plot(ctimes, tempHcrFog, label='HcrFog temp (C)', color='blue', linewidth=1)
    ax2.plot(ctimes, tempTailcone, label='Tailcone temp (C)', color='green', linewidth=1)
    ax2.plot(ctimes, temp30OffsetBy10, label='TempOffsetBy10', color='magenta', linewidth=1)

    ax3.plot(ctimes, ias, \
             label='Indicated Airspeed', color='green', linewidth=1)

    ax4.plot(ctimes, gvMass10000Kg, \
             label='GV mass (*10T)', color='green', linewidth=2)
    ax4.plot(ctimes, aoaSm, \
             label='Angle of attack', color='red', linewidth=1)
    ax4.plot(ctimes, accelNormSm, \
             label='accelNorm', color='orange', linewidth=1)

    configTimeAxis(ax1, -9999, -9999, "Altitude", 'upper center')
    configTimeAxis(ax2, -9999, -9999, "Temp (C)", 'lower center')
    configTimeAxis(ax3, -9999, -9999, "Airspeed (m/s)", 'upper center')
    configTimeAxis(ax4, -9999, -9999, "AOA, G-Accel, Mass", 'upper right')
    # configTimeAxis(ax2, -50, 50, "Temp (C)", 'lower center')
    # configTimeAxis(ax3, 50, 200, "Airspeed (m/s)", 'upper center')
    # configTimeAxis(ax4, 0, 6, "AOA, G-Accel, Mass", 'upper right')

    fig.autofmt_xdate()
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.08, left=0.06, right=0.97, top=0.96)

    fig.suptitle("FLIGHT OVERVIEW - " + str(startTime) + " to " + str(endTime))

    return

########################################################################
# Plot INS and diffs

def doPlotPitchRoll():

    # set up plots

    widthIn = float(options.mainWidthMm) / 25.4
    htIn = float(options.mainHeightMm) / 25.4
    
    global figNum
    fig = plt.figure(figNum, (widthIn, htIn))
    figNum = figNum + 1
    
    ax1 = fig.add_subplot(4,1,1,xmargin=0.0)
    ax2 = fig.add_subplot(4,1,2,xmargin=0.0)
    ax3 = fig.add_subplot(4,1,3,xmargin=0.0)
    ax4 = fig.add_subplot(4,1,4,xmargin=0.0)
    
    #ax1.plot(ctimes, pitch3Sm, label='PitchIns3', color='orange', linewidth=1)
    #ax1.plot(ctimes, pitch2Sm, label='PitchIns2', color='blue', linewidth=1)
    ax1.plot(ctimes, pitchSm, label='PitchIns1', color='green', linewidth=1)
    ax1.plot(ctimes, pitchCmSm, label='PitchHcrFog', color='red', linewidth=1)
    
    #ax2.plot(ctimes, roll3Sm, label='RollIns3', color='orange', linewidth=1)
    #ax2.plot(ctimes, roll2Sm, label='RollIns2', color='green', linewidth=1)
    ax2.plot(ctimes, rollSm, label='RollIns1', color='green', linewidth=1)
    ax2.plot(ctimes, rollCmSm, label='RollHcrFog', color='red', linewidth=1)
    
    #ax3.plot(ctimes, pitchDiffCm2Sm, \
    #         label='pitchDiffCm2Sm', color='magenta', linewidth=1)
    #ax3.plot(ctimes, estPitchDiffSm, \
    #         label='estPitchDiff', color='blue', linewidth=1)
    ax3.plot(ctimes, pitchDiffSm, \
             label='pitchDiffHcrFog', color='blue', linewidth=1)
    #ax3.plot(ctimes, pitchDiff2Sm, \
    #         label='pitchDiffIns2', color='red', linewidth=1)
    #ax3.plot(ctimes, pitchDiff3Sm, \
    # label='pitchDiffIns3', color='black', linewidth=1)


    ax4.plot(ctimes, rollDiffSm, \
             label='rollDiffHcrFog', color='blue', linewidth=1)
    #ax4.plot(ctimes, rollDiff2Sm, \
    #         label='rollDiffIns2', color='red', linewidth=1)
    #ax4.plot(ctimes, rollDiff3Sm, \
    #         label='rollDiffIns3', color='black', linewidth=1)

    #configTimeAxis(ax1, -9999, -9999, "Pitch", 'upper right')
    configTimeAxis(ax2, -9999, -9999, "Roll", 'upper right')
    #configTimeAxis(ax3, -9999, -9999, "PitchDiffs", 'upper right')
    #configTimeAxis(ax4, -9999, -9999, "RollDiffs", 'upper right')

    configTimeAxis(ax1, -2, 4, "Pitch", 'upper right')
    # configTimeAxis(ax2, -1.5, 4, "Roll", 'upper right')
    configTimeAxis(ax3, -1, 3, "PitchDiffs", 'upper right')
    configTimeAxis(ax4, -2, 2, "RollDiffs", 'upper right')

    fig.autofmt_xdate()
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.08, left=0.06, right=0.97, top=0.96)

    # title name

    fig.suptitle("ROLL and PITCH - " + str(startTime) + " to " + str(endTime))

    return

########################################################################
# Plot the diffs

def doPlotDiffs():

    # set up plots

    widthIn = float(options.mainWidthMm) / 25.4
    htIn = float(options.mainHeightMm) / 25.4
    
    global figNum
    fig = plt.figure(figNum, (widthIn, htIn))
    figNum = figNum + 1
    
    ax1 = fig.add_subplot(4,1,1,xmargin=0.0)
    ax2 = fig.add_subplot(4,1,2,xmargin=0.0)
    ax3 = fig.add_subplot(4,1,3,xmargin=0.0)
    ax4 = fig.add_subplot(4,1,4,xmargin=0.0)
    
    ax1.plot(ctimes, driftDiffSm, \
             label='driftDiff', color='black', linewidth=1)
    ax1.plot(ctimes, hdgDiffSm, \
             label='hdgDiff', color='blue', linewidth=1)
    ax1.plot(ctimes, vVelDiffSm, \
             label='vVelDiff', color='green', linewidth=1)
    ax1.plot(ctimes, altDiff, \
             label='altDiff', color='red', linewidth=2)

    ax2.plot(ctimes, estPitchDiffSm, \
             label='estPitchDiff', color='black', linewidth=1)
    # ax2.plot(ctimes, altGps / 1000.0, \
    #          label='Alt(km)', color='gray', linewidth=2)
    #ax2.plot(ctimes, gvMass10000Kg, \
    #         label='Mass (10T)', color='brown', linewidth=1)
    ax2.plot(ctimes, accelNormSm, \
             label='accelNorm', color='orange', linewidth=1)
    ax2.plot(ctimes, aoa2, \
             label='aoa/2', color='red', linewidth=1)
    ax2.plot(ctimes, pitchDiffSm, \
             label='pitchDiff', color='green', linewidth=1)
    ax2.plot(ctimes, rollDiffSm, \
             label='rollDiff', color='blue', linewidth=1)
    #ax2.plot(ctimes, temp30OffsetBy10, label='TempOffsetBy10', color='magenta', linewidth=1)
    

    ax3.plot(ctimes, surfaceVelSm, \
             label='surfaceVel', color='blue', linewidth=1)

    ax4.plot(ctimes, elevErrSm, label='ElevErr', color='red', linewidth=1)
    ax4.plot(ctimes, tiltErrSm, label='TiltErr', color='blue', linewidth=1)

    # configTimeAxis(ax1, -5, 5.0, "diffs", 'upper right')
    # configTimeAxis(ax2, -2, 3, "diffs", 'upper right')
    configTimeAxis(ax1, -9999, -9999, "diffs", 'upper right')
    configTimeAxis(ax2, -9999, -9999, "diffs", 'upper right')
    configTimeAxis(ax3, -10, 10, "SurfaceVel", 'upper right')
    configTimeAxis(ax4, -2, 2, "Err", 'upper right')
    
    # configTimeAxis(ax2, -0.5, 2.5, "diffs", 'upper right')
    # configTimeAxis(ax3, -9999, -9999, "SurfaceVel", 'upper right')
    # configTimeAxis(ax4, -0.25, 0.25, "Err", 'upper right')
    
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.08, left=0.06, right=0.97, top=0.96)

    fig.suptitle("DIFFS GV minus HCR FOG - " + str(startTime) + " to " + str(endTime))

    return

########################################################################
# Produce the 2-D histograms

def doPlot2DHist():

    # set up plots

    widthIn = float(options.mainWidthMm) / 25.4
    htIn = float(options.mainHeightMm) / 25.4
    
    global figNum
    fig = plt.figure(figNum, (widthIn, htIn))
    figNum = figNum + 1
    
    ax1 = fig.add_subplot(2,2,1,xmargin=0.0)
    ax2 = fig.add_subplot(2,2,2,xmargin=0.0)
    ax3 = fig.add_subplot(2,2,3,xmargin=0.0)
    ax4 = fig.add_subplot(2,2,4,xmargin=0.0)

    # pitch offset vs G force

    validAccelNorm = accelNormSm[np.isfinite(accelNormSm)]
    validPitchDiff = pitchDiffSm[np.isfinite(pitchDiffSm)]
    validRollDiff = rollDiffSm[np.isfinite(rollDiffSm)]
    validAoa = aoaSm[np.isfinite(aoaSm)]
    validTemp30Offset = temp30OffsetSm[np.isfinite(temp30OffsetSm)]

    (counts1, xedges1, yedges1, Image1) = ax1.hist2d(validAccelNorm, validPitchDiff, 
                                                     bins=100, norm=mpl.colors.LogNorm())
    cbar1 = fig.colorbar(mappable=Image1, ax=ax1)
    cbar1.set_label("count")
    ax1.set_xlabel('G-acceleration normal')
    ax1.set_ylabel('Pitch offset')

    # roll offset vs G force

    (counts2, xedges2, yedges2, Image2) = ax2.hist2d(validAccelNorm, validRollDiff,
                                                     bins=100, norm=mpl.colors.LogNorm())
    cbar2 = fig.colorbar(mappable=Image2, ax=ax2)
    cbar2.set_label("count")
    ax2.set_xlabel('G-acceleration normal')
    ax2.set_ylabel('Roll offset')

    # pitch offset vs temp diff

    (counts3, xedges3, yedges3, Image3) = ax3.hist2d(validTemp30Offset, validPitchDiff, 
                                                     bins=100, norm=mpl.colors.LogNorm())
    cbar3 = fig.colorbar(mappable=Image3, ax=ax3)
    cbar3.set_label("count")
    ax3.set_xlabel('30 - hcrfog temp')
    ax3.set_ylabel('Pitch offset')

    # roll offset vs Angle of attack

    (counts4, xedges4, yedges4, Image4) = ax4.hist2d(validAoa, validRollDiff, 
                                                     bins=100, norm=mpl.colors.LogNorm())
    cbar4 = fig.colorbar(mappable=Image4, ax=ax4)
    cbar4.set_label("count")
    ax4.set_xlabel('AOA')
    ax4.set_ylabel('Roll offset')

    fig.tight_layout()
    fig.subplots_adjust(bottom=0.08, left=0.06, right=0.97, top=0.96)

    fig.suptitle("2D histograms - " + str(startTime) + " to " + str(endTime))

    return

########################################################################
# Plot angles

def doPlotRadarAngles():

    # set up plots

    widthIn = float(options.mainWidthMm) / 25.4
    htIn = float(options.mainHeightMm) / 25.4
    
    global figNum
    fig = plt.figure(figNum, (widthIn, htIn))
    figNum = figNum + 1
    
    ax1 = fig.add_subplot(4,1,1,xmargin=0.0)
    ax2 = fig.add_subplot(4,1,2,xmargin=0.0)
    ax3 = fig.add_subplot(4,1,3,xmargin=0.0)
    ax4 = fig.add_subplot(4,1,4,xmargin=0.0)
    
    ax1.plot(ctimes, azimuth, label='Azimuth', color='red', linewidth=1)

    ax2.plot(ctimes, elevErrSm, label='ElevErr', color='blue', linewidth=1)

    ax3.plot(ctimes, rotation, label='Rotation', color='green', linewidth=1)

    ax4.plot(ctimes, tilt, label='Tilt', color='black', linewidth=1)

    configTimeAxis(ax1, -9999, -9999, "Azimuth", 'upper right')
    configTimeAxis(ax2, -0.25, 0.25, "Elevation error", 'upper right')
    configTimeAxis(ax3, -9999, -9999, "Rotation", 'upper right')
    configTimeAxis(ax4, -9999, -9999, "Tilt", 'upper right')

    fig.autofmt_xdate()
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.08, left=0.06, right=0.97, top=0.96)

    # title name

    fig.suptitle("RADAR ANGLES - " + str(startTime) + " to " + str(endTime))

    return

########################################################################
# Plot the estimated pitch difference, using surface vel to estimate it

def doPlotEstPitchDiff():

    # set up plots

    widthIn = float(options.mainWidthMm) / 25.4
    htIn = float(options.mainHeightMm) / 25.4
    
    global figNum
    fig = plt.figure(figNum, (widthIn, htIn))
    figNum = figNum + 1
    
    ax1 = fig.add_subplot(2,1,1,xmargin=0.0)
    ax2 = fig.add_subplot(2,1,2,xmargin=0.0)
    
    ax1.plot(ctimes, estPitchDiffSm, \
             label='estPitchDiff', color='red', linewidth=1)
    ax1.plot(ctimes, pitchDiffSm, \
             label='pitchDiff', color='green', linewidth=1)
    ax1.plot(ctimes, rollDiffSm, \
             label='rollDiff', color='blue', linewidth=1)
    
    ax2.plot(ctimes, surfaceVelSm, \
             label='surfaceVel', color='blue', linewidth=1)

    configTimeAxis(ax1, -0.5, 1, "diffs", 'upper right')
    configTimeAxis(ax2, -1.5, 1, "SurfaceVel", 'upper right')
    
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.08, left=0.06, right=0.97, top=0.96)

    fig.suptitle("ESTIMATED PITCH DIFF: GV-HCRFOG - " + str(startTime) + " to " + str(endTime))
    
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
# Create the comparison file, running AcGeorefCompare

def createCompFile():

    projDir = os.environ['PROJ_DIR']
    paramsPath = projDir + '/qc/params/AcGeorefCompare.iwg1_vs_10hz'
    startTimeStr = '"' + \
                   str(startTime.year) + ' ' + str(startTime.month) + ' ' + \
                   str(startTime.day) + ' ' + str(startTime.hour) + ' ' + \
                   str(startTime.minute) + ' ' + str(startTime.second) + '"'
    endTimeStr = '"' + \
                 str(endTime.year) + ' ' + str(endTime.month) + ' ' + \
                 str(endTime.day) + ' ' + str(endTime.hour) + ' ' + \
                 str(endTime.minute) + ' ' + str(endTime.second) + '"'

    command = 'AcGeorefCompare -params ' + paramsPath + \
              ' -start ' + startTimeStr + ' -end ' + endTimeStr + \
              ' > ' + compFilePath

    runCommand(command)


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

