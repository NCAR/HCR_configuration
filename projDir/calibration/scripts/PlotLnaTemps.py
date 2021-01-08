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
from numpy import linalg, array, ones
import scipy
from scipy import odr
import matplotlib.pyplot as plt
from matplotlib import dates
import math
import datetime

def main():

#   globals

    global options
    global debug
    global startTime
    global endTime

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
                      default=250,
                      help='Width of figure in mm')
    parser.add_option('--height',
                      dest='figHeightMm',
                      default=140,
                      help='Height of figure in mm')
    parser.add_option('--start',
                      dest='startTime',
                      default='2000 01 01 00 00 00',
                      help='Start time for plots')
    parser.add_option('--end',
                      dest='endTime',
                      default='2100 01 01 00 00 00',
                      help='End time for plots')
    parser.add_option('--pwrMin',
                      dest='pwrMinDbm',
                      default='-9999.0',
                      help='Minimum power on plots (dBm)')
    parser.add_option('--pwrMax',
                      dest='pwrMaxDbm',
                      default='-9999.0',
                      help='Minimum power on plots (dBm)')
    parser.add_option('--lenMean',
                      dest='lenMean',
                      default=1,
                      help='Len of moving mean filter')
    parser.add_option('--lagSecs',
                      dest='lagSecs',
                      default=0,
                      help='Correct for lag between temp and power measurements')
    parser.add_option('--hOnly',
                      dest='hOnly', default=False,
                      action="store_true",
                      help='Only plot H channel data sets')
    parser.add_option('--vOnly',
                      dest='vOnly', default=False,
                      action="store_true",
                      help='Only plot V channel data sets')
    parser.add_option('--lag1',
                      dest='lag1', default=False,
                      action="store_true",
                      help='Plot lag1 coherent power instead of lag0 power')
    
    (options, args) = parser.parse_args()
    
    if (options.verbose):
        options.debug = True

    year, month, day, hour, minute, sec = options.startTime.split()
    startTime = datetime.datetime(int(year), int(month), int(day),
                                  int(hour), int(minute), int(sec))

    year, month, day, hour, minute, sec = options.endTime.split()
    endTime = datetime.datetime(int(year), int(month), int(day),
                                int(hour), int(minute), int(sec))
    if (options.debug):
        print("Running %prog", file=sys.stderr)
        print("  filePath: ", options.filePath, file=sys.stderr)
        print("  startTime: ", startTime, file=sys.stderr)
        print("  endTime: ", endTime, file=sys.stderr)
        print("  lenMean: ", options.lenMean, file=sys.stderr)
        print("  lagSecs: ", options.lagSecs, file=sys.stderr)
        print("  hOnly: ", options.hOnly, file=sys.stderr)
        print("  vOnly: ", options.vOnly, file=sys.stderr)
        print("  lag1: ", options.lag1, file=sys.stderr)

    if (options.hOnly and options.vOnly):
        print("ERROR - cannot set both --hOnly and --vOnly", file=sys.stderr)
        sys.exit(1)

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
                if (thisTime < startTime or thisTime > endTime):
                    break
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

    for ii in range(0, int(window / 2) + 1):
        sma[ii] = values[ii]
        sma[len(sma)-ii-1] = values[len(sma)-ii-1]

    return sma

########################################################################
# define funtion for linear fit

def flinear(B, x):
    '''Linear function y = m*x + b'''
    # B is a vector of the parameters.
    # x is an array of the current x values.
    # x is in the same format as the x passed to Data or RealData.
    #
    # Return an array in the same format as y passed to Data or RealData.
    return B[0]*x + B[1]

########################################################################
# Plot

def doPlot(colHdrs, obsTimes, colData):

    lenMeanFilter = int(options.lenMean)

    fileName = os.path.basename(options.filePath)
    titleStr = "File: " + fileName
    hfmt = dates.DateFormatter('%y/%m/%d')

    # times
    
    obstimes = np.array(obsTimes).astype(datetime.datetime)

    # power
    

    powerHc = np.array(colData["Hc"]).astype(np.double)
    powerVc = np.array(colData["Vc"]).astype(np.double)
    lag1Hc = np.array(colData["Lag1Hc"]).astype(np.double)
    lag1Vc = np.array(colData["Lag1Vc"]).astype(np.double)

    if (options.hOnly):
        minPwr = min(powerHc)
        maxPwr = max(powerHc)
        if (options.lag1):
            minPwr = min(minPwr, min(lag1Hc))
    elif (options.vOnly):
        minPwr = min(powerVc)
        maxPwr = max(powerVc)
        if (options.lag1):
            minPwr = min(minPwr, min(lag1Vc))
    else:
        minPwr = min(min(powerHc), min(powerVc))
        maxPwr = max(max(powerHc), max(powerVc))
        if (options.lag1):
            minPwr = min(minPwr, min(lag1Hc))
            minPwr = min(minPwr, min(lag1Vc))

    rangePwr = maxPwr - minPwr

    # temps with moving average

    tempH = np.array(colData["HLnaTemp"]).astype(np.double)
    tempV = np.array(colData["VLnaTemp"]).astype(np.double)

    lagSecs = int(options.lagSecs)
    if (lagSecs > 0):
        for ii in range(0, len(tempH) - lagSecs):
            tempH[ii] = tempH[ii + lagSecs]
            tempV[ii] = tempV[ii + lagSecs]

    tempH = movingAverage(tempH, lenMeanFilter)
    tempV = movingAverage(tempV, lenMeanFilter)

    if (options.hOnly):
        minTemp = min(tempH)
        maxTemp = max(tempH)
    elif (options.vOnly):
        minTemp = min(tempV)
        maxTemp = max(tempV)
    else:
        minTemp = min(min(tempH), min(tempV))
        maxTemp = max(max(tempH), max(tempV))

    rangeTemp = maxTemp - minTemp

    # set up plot structure

    widthIn = float(options.figWidthMm) / 25.4
    htIn = float(options.figHeightMm) / 25.4

    fig1 = plt.figure(1, (widthIn, htIn))
    ax1 = fig1.add_subplot(2,1,1,xmargin=0.0)
    ax2 = fig1.add_subplot(2,1,2,xmargin=0.0)
    ax1.set_xlim([obstimes[0], obstimes[-1]])
    ax2.set_xlim([obstimes[0], obstimes[-1]])

    # axis 1 - power
    
    if (not options.vOnly):
        ax1.plot(obstimes, powerHc, \
                 label = 'powerHc', linewidth=1, color='blue')
        if (options.lag1):
            ax1.plot(obstimes, lag1Hc, \
                     label = 'lag1Hc', linewidth=2, linestyle=':', color='orange')
        
    if (not options.hOnly):
        ax1.plot(obstimes, powerVc, \
                 label = 'powerVc', linewidth=1, color='red')
        if (options.lag1):
            ax1.plot(obstimes, lag1Vc, \
                     label = 'lag1Vc', linewidth=2, linestyle=':', color='green')

    # axis 2 - temps
    
    if (not options.vOnly):
        ax2.plot(obstimes, tempH, \
                 label = 'tempH', linewidth=1, color='blue')

    if (not options.hOnly):
        ax2.plot(obstimes, tempV, \
                 label = 'tempV', linewidth=1, color='red')

    # labels

    ax1.set_title("Received power (dBm)", fontsize=10)
    ax2.set_title("LNA Temp (C)", fontsize=10)
    ax1.set_ylim(minPwr - rangePwr * 0.1, maxPwr + rangePwr * 0.3)
    ax2.set_ylim(minTemp - rangeTemp * 0.1, maxTemp + rangeTemp * 0.2)

    configureTimeAxis(ax1, float(options.pwrMinDbm), float(options.pwrMaxDbm),
                      "Power", 'upper left')
    configureTimeAxis(ax2, -9999.0, -9999.0,
                      "Temps", 'upper left')

    title1 = "HCR LNA Power/Temp Time Series " + \
             startTime.isoformat(sep=" ") + \
             " to " + endTime.isoformat(sep=" ")
    fig1.suptitle(title1, fontsize=12)
    fig1.autofmt_xdate()
    fig1.subplots_adjust(bottom=0.12, left=0.10, right=0.95, top=0.90)

    # save ax1/2 plot to file

    if (options.hOnly):
        saveLabel = '.hOnly'
    elif (options.vOnly):
        saveLabel = '.vOnly'
    else:
        saveLabel = ''

    homeDir = os.environ['HOME']
    saveDir = os.path.join(homeDir, 'Downloads')
    saveDir = os.path.join(saveDir, 'images')
    saveName1 = 'power_vs_temp_timeseries.' + \
               startTime.isoformat() + "-" + endTime.isoformat() + \
               saveLabel + '.png'
    savePath1 = os.path.join(saveDir, saveName1)
    print("  saving ax1/2 figure to path: ", savePath1, file=sys.stderr)
    plt.savefig(savePath1, pad_inches=0.0)

    # Plot of temp vs gain, with orthogonal linear fits

    # Horiz

    fig3 = plt.figure(2, (int(widthIn/1.5), int(htIn/1.5)))
    ax3 = fig3.add_subplot(1,1,1,xmargin=1.0, ymargin=1.0)

    if (not options.vOnly):
        ax3.plot(tempH, powerHc, ".", color = 'blue')
        if (options.lag1):
            ax3.plot(tempH, lag1Hc, "x", color = 'orange')

    linear = odr.Model(flinear)
    dataH = odr.Data(tempH, powerHc)
    odrH = odr.ODR(dataH, linear, beta0=[1.0, 2.0])
    outputH = odrH.run()
    #outputH.pprint()
    slopeH = outputH.beta[0]
    interceptH = outputH.beta[1]

    #AH = array([tempH, ones(len(tempH))])
    # obtaining the fit, ww[0] is slope, ww[1] is intercept
    #wwH = linalg.lstsq(AH.T, powerHc)[0]

    regrXH = []
    regrYH = []
    minTempH = min(tempH) - 5.0
    maxTempH = max(tempH) + 5.0
    minPwrH = slopeH * minTempH + interceptH
    maxPwrH = slopeH * maxTempH + interceptH
    regrXH.append(minTempH)
    regrXH.append(maxTempH)
    regrYH.append(minPwrH)
    regrYH.append(maxPwrH)
    labelH = "Gain slope H = " + ("%.3f" % slopeH)
    if (not options.vOnly):
        ax3.plot(regrXH, regrYH, linewidth=1, color = 'blue', label=labelH)

    # Vert

    if (not options.hOnly):
        ax3.plot(tempV, powerVc, ".", color = 'red')
        if (options.lag1):
            ax3.plot(tempV, lag1Vc, "x", color = 'green')

    dataV = odr.Data(tempV, powerVc)
    odrV = odr.ODR(dataV, linear, beta0=[1.0, 2.0])
    outputV = odrV.run()
    #outputV.pprint()
    slopeV = outputV.beta[0]
    interceptV = outputV.beta[1]

    #AV = array([tempV, ones(len(tempV))])
    # obtaining the fit, ww[0] is slope, ww[1] is intercept
    #wwV = linalg.lstsq(AV.T, powerVc)[0]

    regrXV = []
    regrYV = []
    minTempV = min(tempV) - 5.0
    maxTempV = max(tempV) + 5.0
    minPwrV = slopeV * minTempV + interceptV
    maxPwrV = slopeV * maxTempV + interceptV
    regrXV.append(minTempV)
    regrXV.append(maxTempV)
    regrYV.append(minPwrV)
    regrYV.append(maxPwrV)
    labelV = "Gain slope V = " + ("%.3f" % slopeV)
    if (not options.hOnly):
        ax3.plot(regrXV, regrYV, linewidth=1, color = 'red', label=labelV)

    minTemp2 = min(min(tempH), min(tempV))
    maxTemp2 = max(max(tempH), max(tempV))
    rangeTemp2 = maxTemp2 - minTemp2

    ax3.set_xlim(minTemp2 - rangeTemp2 * 0.2, maxTemp2 + rangeTemp2 * 0.2)
    ax3.set_ylim(minPwr - rangeTemp2 * 0.2, maxPwr + rangeTemp2 * 0.25)

    title3 = "Power vs Temp " + \
             startTime.isoformat(sep=" ") + \
             " to " + endTime.isoformat(sep=" ")
    ax3.set_title(title3, fontsize=10)
    labelX = "Temp(C) lagged " + ("%d" % lagSecs) + " secs"
    ax3.set_xlabel(labelX)
    ax3.set_ylabel("Measured power(dBm)")
    
    legend3 = ax3.legend(loc="upper left", ncol=2)
    fig3.subplots_adjust(bottom=0.15, left=0.13, right=0.95, top=0.90)

    # save ax3 plot to file

    saveName3 = 'power_vs_temp_xyplot.' + \
               startTime.isoformat() + "-" + endTime.isoformat() + \
               saveLabel + '.png'
    savePath3 = os.path.join(saveDir, saveName3)
    print("  saving ax3 figure to path: ", savePath3, file=sys.stderr)
    plt.savefig(savePath3, pad_inches=0.0)

    # show

    plt.tight_layout()
    plt.show()


########################################################################
# initialize legends etc

def configureTimeAxis(ax, miny, maxy, ylabel, legendLoc):
    
    legend = ax.legend(loc=legendLoc, ncol=6)
    for label in legend.get_texts():
        label.set_fontsize('x-small')
    ax.set_xlabel("Time")
    ax.set_ylabel(ylabel)
    ax.grid(True)
    if (miny > -9990 and maxy > -9990):
        ax.set_ylim([miny, maxy])
    hfmt = dates.DateFormatter('%H:%M:%S')
    #ax.xaxis.set_major_locator(dates.HourLocator())
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

