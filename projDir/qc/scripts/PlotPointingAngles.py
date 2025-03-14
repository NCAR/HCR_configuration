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
                      default=4.1,
                      help='Maximum pitch angle (deg)')
    parser.add_option('--pitchDelta',
                      dest='pitchDelta',
                      default=1.0,
                      help='Delta pitch angle in matrix (deg)')
    parser.add_option('--rollMin',
                      dest='rollMin',
                      default=-30.0,
                      help='Minimum roll angle (deg)')
    parser.add_option('--rollMax',
                      dest='rollMax',
                      default=30.1,
                      help='Maximum roll angle (deg)')
    parser.add_option('--rollDelta',
                      dest='rollDelta',
                      default=10.0,
                      help='Delta roll angle in matrix (deg)')
    
    (options, args) = parser.parse_args()

    global _tilt, _rot, _drift
    global _pitchMin, _pitchMax, _pitchDelta
    global _rollMin, _rollMax, _rollDelta

    _tilt = float(options.tilt)
    _rot = float(options.rot)
    _drift = float(options.drift)

    _pitchMin = float(options.pitchMin)
    _pitchMax = float(options.pitchMax)
    _pitchDelta = float(options.pitchDelta)

    _rollMin = float(options.rollMin)
    _rollMax = float(options.rollMax)
    _rollDelta = float(options.rollDelta)

    print("#  tilt: ", _tilt, file=sys.stderr)
    print("#  rot: ", _rot, file=sys.stderr)
    print("#  drift: ", _drift, file=sys.stderr)
    print("#  pitchMin: ", _pitchMin, file=sys.stderr)
    print("#  pitchMax: ", _pitchMax, file=sys.stderr)
    print("#  pitchDelta: ", _pitchDelta, file=sys.stderr)
    print("#  rollMin: ", _rollMin, file=sys.stderr)
    print("#  rollMax: ", _rollMax, file=sys.stderr)
    print("#  rollDelta: ", _rollDelta, file=sys.stderr)

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

    print("# ",
          '{:>10} '.format('pitch'),
          '{:>10} '.format('roll'),
          '{:>10} '.format('reflRot'),
          '{:>10} '.format('reflTilt'),
          file=sys.stderr)

    for pitch in np.arange(_pitchMin, _pitchMax, _pitchDelta):

        for roll in np.arange(_rollMin, _rollMax, _rollDelta):

            reflRot, reflTilt = computeReflAngles(pitch, roll)

            print("  ",
                  '{:10.2f} '.format(pitch),
                  '{:10.2f} '.format(roll),
                  '{:10.2f} '.format(reflRot),
                  '{:10.2f} '.format(reflTilt),
                  file=sys.stderr)

########################################################################
# compute reflector angles
# Returns rotation and tilt for the reflector
# This is from MotionControl.cpp.

def computeReflAngles(pitch, roll):

    sinPitch = math.sin(math.radians(pitch))
    cosPitch = math.cos(math.radians(pitch))

    sinRoll = math.sin(math.radians(roll))
    cosRoll = math.cos(math.radians(roll))

    sinDrift = math.sin(math.radians(_drift))
    cosDrift = math.cos(math.radians(_drift))

    # Track relative coordinates - desired beam position

    sinRot = math.sin(math.radians(_rot))
    cosRot = math.cos(math.radians(_rot))

    sinTilt = math.sin(math.radians(_tilt))
    cosTilt = math.cos(math.radians(_tilt))

    # Convert to track relative Cartesian coordinates

    x_t = cosTilt * sinRot
    y_t = sinTilt
    z_t = cosTilt * cosRot

    #print("#  x_t, y_t, z_t: ", x_t, y_t, z_t, file=sys.stderr)

    # Convert to pod relative Cartesian coordinates - adjusted beam position

    x_a = \
          x_t * (cosDrift * cosRoll - sinDrift * sinPitch * sinRoll) + \
          y_t * (sinDrift + cosDrift * sinPitch * sinRoll) + \
          -z_t * cosPitch * sinRoll

    y_a = \
          -x_t * sinDrift * cosPitch + \
          y_t * cosDrift * cosPitch + \
          z_t * sinPitch

    z_a = \
          x_t * (cosDrift * sinRoll + sinDrift * sinPitch * sinRoll) + \
          y_t * (sinDrift * sinRoll - cosDrift * sinPitch * cosRoll) + \
          z_t * cosPitch * cosRoll

    #print("#  x_a, y_a, z_a: ", x_a, y_a, z_a, file=sys.stderr)

    # Convert from pod relative Cartesian coordinates to polar coordinates
    # and save the adjusted rotation and tilt angles.

    reflTilt = math.degrees(math.asin(y_a))

    # KLUGE: The algorithm above isn't really right. As compensation for now,
    # just change the sign of the corrected tilt if the desired rotation angle
    # rots downward, i.e., if its cosine is less than zero.
    if (cosRot < 0.0):
        reflTilt = reflTilt * -1

    # rotation needs to be in the range of 0-360

    reflRot = math.degrees(math.atan2(x_a, z_a))
    #if (reflRot < 0):
    #    reflRot = reflRot + 360.0

    return (reflRot, reflTilt)

########################################################################
# compute (elevation, azimuth) from rotation and tilt, attitude
# For a Y-prime radar e.g. HCR

def computeAzElYPrime(pitch, roll, hdg, rot, tilt):

    # precompute sin/cos
    
    sinPitch = math.sin(math.radians(pitch))
    cosPitch = math.cos(math.radians(pitch))

    sinRoll = math.sin(math.radians(roll))
    cosRoll = math.cos(math.radians(roll))

    sinHdg = math.sin(math.radians(hdg))
    cosHdg = math.cos(math.radians(hdg))

    sinRot = math.sin(math.radians(rot))
    cosRot = math.cos(math.radians(rot))
    
    sinTilt = math.sin(math.radians(tilt))
    cosTilt = math.cos(math.radians(tilt))

    # compute unit vector relative to aircraft

    x_a = sinRot * cosTilt
    y_a = sinTilt
    z_a = cosRot * cosTilt

    # compute matrix elements after multiplication
    # for 3 axis transformation

    m11 = cosHdg * cosRoll + sinHdg * sinPitch * sinRoll
    m12 = sinHdg * cosPitch
    m13 = cosHdg * sinRoll - sinHdg * sinPitch * cosRoll

    m21 = -sinHdg * cosRoll + cosHdg * sinPitch * sinRoll
    m22 = cosHdg * cosPitch
    m23 = -sinHdg * sinRoll - cosHdg * sinPitch * cosRoll

    m31 = -cosPitch * sinRoll
    m32 = sinPitch
    m33 = cosPitch * cosRoll

    # Compute unit vector in earth coords

    xx = m11 * x_a + m12 * y_a * m13 * z_a
    yy = m21 * x_a + m22 * y_a * m23 * z_a
    zz = m31 * x_a + m32 * y_a * m33 * z_a

    # compute az and el

    az = math.degrees(math.atan2(xx, yy))
    el = math.degrees(math.asin(zz))

    print("#  pitch: ", pitch, file=sys.stderr)
    print("#  roll: ", roll, file=sys.stderr)
    print("#  hdg: ", hdg, file=sys.stderr)
    print("#  rot: ", rot, file=sys.stderr)
    print("#  tilt: ", tilt, file=sys.stderr)

    print("#  x_a: ", x_a, file=sys.stderr)
    print("#  y_a: ", y_a, file=sys.stderr)
    print("#  z_a: ", z_a, file=sys.stderr)

    print("#  xx: ", x_a, file=sys.stderr)
    print("#  yy: ", yy, file=sys.stderr)
    print("#  zz: ", zz, file=sys.stderr)

    print("#  az: ", az, file=sys.stderr)
    print("#  el: ", el, file=sys.stderr)

    return (el, az)

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

