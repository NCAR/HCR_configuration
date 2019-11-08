#!/usr/bin/env python

#===========================================================================
#
# Produce plots for OTREC RF09 data
#
#===========================================================================

from __future__ import print_function
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

def main():

    # globals

    global options
    global debug
    global figNum
    figNum = 0

    global appName
    appName = os.path.basename(__file__)

    # parse the command line

    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option('--debug',
                      dest='debug', default=False,
                      action="store_true",
                      help='Set debugging on')
    parser.add_option('--tilt',
                      dest='tilt',
                      default=0.0,
                      help='Desired tilt angle (deg)')
    parser.add_option('--rot',
                      dest='rot',
                      default=-90.0,
                      help='Desired rotation angle (deg)')
    parser.add_option('--drift',
                      dest='drift',
                      default=0.0,
                      help='Desired drift angle (deg)')
    parser.add_option('--pitchMin',
                      dest='pitchMin',
                      default=-2.0,
                      help='Minimum pitch angle (deg)')
    parser.add_option('--pitchMax',
                      dest='pitchMax',
                      default=4.0,
                      help='Maximum pitch angle (deg)')
    parser.add_option('--pitchDelta',
                      dest='pitchDelta',
                      default=0.5,
                      help='Delta pitch angle in matrix (deg)')
    parser.add_option('--rollMin',
                      dest='rollMin',
                      default=-30.0,
                      help='Minimum roll angle (deg)')
    parser.add_option('--rollMax',
                      dest='rollMax',
                      default=30.0,
                      help='Maximum roll angle (deg)')
    parser.add_option('--rollDelta',
                      dest='rollDelta',
                      default=5.0,
                      help='Delta roll angle in matrix (deg)')
    
    (options, args) = parser.parse_args()
    
    if (options.debug == True):
        print("  tilt: ", options.tilt, file=sys.stderr)
        print("  rot: ", options.rot, file=sys.stderr)
        print("  drift: ", options.drift, file=sys.stderr)
        print("  pitchMin: ", options.pitchMin, file=sys.stderr)
        print("  pitchMax: ", options.pitchMax, file=sys.stderr)
        print("  pitchDelta: ", options.pitchDelta, file=sys.stderr)
        print("  rollMin: ", options.rollMin, file=sys.stderr)
        print("  rollMax: ", options.rollMax, file=sys.stderr)
        print("  rollDelta: ", options.rollDelta, file=sys.stderr)

    # compute angles, and write them out

    computeAngles()

    # render the plots

    #doPlotPitchRollHeading()
    #doPlotVelocity()
    #doPlotPosition()

    # show them

    #plt.show()

    sys.exit(0)
    
########################################################################
# Compute the angles, and write them out

def computeAngles():

    for pitch in np.arange(options.pitchMin, options.pitchMax, options.pitchDelta):

        for roll in np.arange(options.rollMin, options.rollMax, options.rollDelta):

            print(" pitch, roll: ", pitch, ", ", roll, file=sys.stderr)


########################################################################
# Set up arrays for plotting

def loadDataArrays(compData, compTimes):

    filtLen = int(options.filtLen)
    
    # set up arrays

    global ctimes

    ctimes = np.array(compTimes).astype(datetime.datetime)

    global lat, lon, alt
    lat = movingAverage(np.array(compData["Lat"]).astype(np.double), filtLen)
    lon = movingAverage(np.array(compData["Lon"]).astype(np.double), filtLen)
    alt = movingAverage(np.array(compData["Alt"]).astype(np.double), filtLen)

    global velNorth, velEast, velUp
    velNorth = movingAverage(np.array(compData["VelNorth"]).astype(np.double), filtLen)
    velEast = movingAverage(np.array(compData["VelEast"]).astype(np.double), filtLen)
    velUp = movingAverage(np.array(compData["VelUp"]).astype(np.double), filtLen)

    global pitch, roll, heading
    pitch = movingAverage(np.array(compData["Pitch"]).astype(np.double), filtLen)
    roll = movingAverage(np.array(compData["Roll"]).astype(np.double), filtLen)
    heading = movingAverage(np.array(compData["Heading"]).astype(np.double), filtLen)

    global gndSpdKnots, vertSpdFps
    #for (vn in velNorth, ve in velEast):
    #    print >>sys.stderr, "vn, ve: ", vn, ve

    gndSpdKnots = math.sqrt((velNorth * velNorth) + math.sqrt(velEast * velEast))
    
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

    fig.suptitle("HCR PITCH/ROLL/HDG - " + str(startTime) + " to " + str(endTime))

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

    fig.suptitle("HCR VELOCITIES - " + str(startTime) + " to " + str(endTime))

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
    configTimeAxis(ax1, 39, 41, "latitude", 'upper right')
    configTimeAxis(ax2, -106, -100, "longitude", 'upper right')
    configTimeAxis(ax3, -9999, -9999, "altitude", 'upper right')

    fig.autofmt_xdate()
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.08, left=0.06, right=0.97, top=0.96)
    
    # title name

    fig.suptitle("HCR LAT/LON/ALT POSITION - " + str(startTime) + " to " + str(endTime))

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

