#!/usr/bin/env python

#===========================================================================
#
# Produce plots for LNA temp analysis
#
#===========================================================================

from __future__ import print_function

import os
import sys
import subprocess
from optparse import OptionParser
import numpy as np
import numpy.ma as ma
from numpy import convolve
import matplotlib.pyplot as plt
from matplotlib import dates
import math
import datetime

def main():

#   globals

    global options
    global debug

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
                      dest='filePath',
                      default='../data/lna_temp/TsPrint.lna_temp.20201105.txt',
                      help='TsPrint output file containing powers and LNA temps')
    parser.add_option('--title',
                      dest='title',
                      default='HCR LNA Temperature analysis',
                      help='Title for plot')
    parser.add_option('--width',
                      dest='figWidthMm',
                      default=400,
                      help='Width of figure in mm')
    parser.add_option('--height',
                      dest='figHeightMm',
                      default=320,
                      help='Height of figure in mm')
    
    (options, args) = parser.parse_args()
    
    if (options.verbose):
        options.debug = True

    if (options.debug):
        print("Running %prog", file=sys.stderr)
        print("  filePath: ", options.filePath, file=sys.stderr)

    # read in column headers for TsPrint output

    iret, colHdrs, colData = readColumnHeaders(options.filePath)
    if (iret != 0):
        sys.exit(-1)

    # read in data

    (obsTimes, colData) = readInputData(options.filePath, colHdrs, colData)

    # render the plot
    
    doPlot(colHdrs, obsTimes, colData)

    sys.exit(0)
    
########################################################################
# Read columm headers for the data
# this is in the first line

def readColumnHeaders(filePath):

    colHeaders = []
    colData = {}

    # read in all lines
    
    fp = open(filePath, 'r')
    lines = fp.readlines()
    fp.close()
    
    # check a line at a time, looking for a comment line
    # with 'time' in it.

    for line in lines:
        
        if (line.find("#") != 0):
            continue

        if (line.find("time") <= 0):
            continue
            
        # header
        colHeaders = line.lstrip("# ").rstrip("\n").split()
        if (options.debug):
            print("colHeaders: ", colHeaders, file=sys.stderr)
    
        for index, var in enumerate(colHeaders, start=0):
            colData[var] = []

        break
        
    return 0, colHeaders, colData

########################################################################
# Read in the data

def readInputData(filePath, colHeaders, colData):

    obsTimes = []

    # open file

    fp = open(filePath, 'r')
    lines = fp.readlines()
    fp.close()

    # process a line at a time, set colData

    for line in lines:

        if (options.verbose):
            print("line==>>", line, "<<==", file=sys.stderr)
            
        if (line.find("#") >= 0):
            continue

        if (line.find("/") < 0):
            continue
            
        # data
        
        data = line.strip().split()

        for index, var in enumerate(colHeaders, start=0):
            if (var == 'time'):
                # comes in as yyyy/mm/dd_hh:mm:ss.fraction
                dateTimeStr = data[index]
                thisTime = decodeDateTime(dateTimeStr)
                obsTimes.append(thisTime)
            else:
                if (isNumber(data[index])):
                    colData[var].append(float(data[index]))
                else:
                    colData[var].append(data[index])

    return obsTimes, colData

########################################################################
# decode date and time

def decodeDateTime(dateTimeStr):

    dateTimeParts = dateTimeStr.split('_')
    dateStr = dateTimeParts[0]
    timeStr = dateTimeParts[1]

    dateParts = dateStr.split("/")
    timeParts = timeStr.split(":")

    year = int(dateParts[0])
    month = int(dateParts[1])
    day = int(dateParts[2])

    hour = int(timeParts[0])
    minute = int(timeParts[1])
    secs = float(timeParts[2])
    sec = int(secs)
    usec = int((secs - sec) * 1000000.0)

    thisTime = datetime.datetime(year, month, day,
                                 hour, minute, sec, usec)
    
    return thisTime

########################################################################
# Check is a number

def isNumber(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

########################################################################
# Moving average filter

def movingAverage(values, window):

    weights = np.repeat(1.0, window)/window
    sma = np.convolve(values, weights, 'same')
    return sma

########################################################################
# Plot

def doPlot(colHdrs, obsTimes, colData):

    fileName = options.filePath
    titleStr = "File: " + fileName
    hfmt = dates.DateFormatter('%y/%m/%d')

    # times
    
    obstimes = np.array(obsTimes).astype(datetime.datetime)

    # power
    
    powerHc = np.array(colData["Hc"]).astype(np.double)
    powerVc = np.array(colData["Vc"]).astype(np.double)

    # temps

    tempH = np.array(colData["HLnaTemp"]).astype(np.double)
    tempV = np.array(colData["VLnaTemp"]).astype(np.double)

    # set up plot structure

    widthIn = float(options.figWidthMm) / 25.4
    htIn = float(options.figHeightMm) / 25.4

    fig1 = plt.figure(1, (widthIn, htIn))
    ax1 = fig1.add_subplot(2,1,1,xmargin=0.0)
    ax2 = fig1.add_subplot(2,1,2,xmargin=0.0)
    #ax3 = fig1.add_subplot(4,1,3,xmargin=0.0)
    #ax4 = fig1.add_subplot(4,1,4,xmargin=0.0)

    ax1.set_xlim([obstimes[0], obstimes[-1]])
    ax2.set_xlim([obstimes[0], obstimes[-1]])
    #ax3.set_xlim([obstimes[0], obstimes[-1]])
    #ax4.set_xlim([obstimes[0], obstimes[-1]])

    # axis 1 - power
    
    ax1.plot(obstimes, powerHc, \
             label = 'powerHc', linewidth=1, color='blue')

    ax1.plot(obstimes, powerVc, \
             label = 'powerVc', linewidth=1, color='red')

    # axis 2 - temps
    
    ax2.plot(obstimes, tempH, \
             label = 'tempH', linewidth=1, color='blue')

    ax2.plot(obstimes, tempV, \
             label = 'tempV', linewidth=1, color='red')

    # legends labels etc

    ax1.set_title("Received power (dBm)", fontsize=12)
    ax2.set_title("LNA Temp (C)", fontsize=12)
    
    configureAxis(ax1, -9999.0, -9999.0, "Power", 'upper left')
    configureAxis(ax2, -9999.0, -9999.0, "Temps", 'upper left')

    fig1.suptitle("HCR LNA Temperature Dependency", fontsize=16)
    fig1.autofmt_xdate()

    plt.tight_layout()
    fig1.subplots_adjust(bottom=0.10, left=0.06, right=0.97, top=0.90)
    plt.show()

########################################################################
# initialize legends etc

def configureAxis(ax, miny, maxy, ylabel, legendLoc):
    
    legend = ax.legend(loc=legendLoc, ncol=6)
    for label in legend.get_texts():
        label.set_fontsize('x-small')
    ax.set_xlabel("Date")
    ax.set_ylabel(ylabel)
    ax.grid(True)
    if (miny > -9990 and maxy > -9990):
        ax.set_ylim([miny, maxy])
    hfmt = dates.DateFormatter('%y/%m/%d')
    ax.xaxis.set_major_locator(dates.HourLocator())
    ax.xaxis.set_major_formatter(hfmt)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(8) 

########################################################################
# Run a command in a shell, wait for it to complete

def runCommand(cmd):

    if (options.debug):
        print("running cmd:",cmd, file=sys.stderr)
    
    try:
        retcode = subprocess.call(cmd, shell=True)
        if retcode < 0:
            print("Child was terminated by signal: ", -retcode, file=sys.stderr)
        else:
            if (options.debug):
                print("Child returned code: ", retcode, file=sys.stderr)
    except OSError as e:
        print("Execution failed:", e, file=sys.stderr)

########################################################################
# Run - entry point

if __name__ == "__main__":
   main()

