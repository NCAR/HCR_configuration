#!/usr/bin/env python

#===========================================================================
#
# Compute antenna angles for mobile platform
# Initially for Yprime radar - e.g. HCR or Eldora
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
    parser.add_option('--az',
                      dest='az',
                      default=0.0,
                      help='Desired azimuth angle (deg)')
    parser.add_option('--el',
                      dest='el',
                      default=-90.0,
                      help='Desired elevation angle (deg)')
    parser.add_option('--tilt',
                      dest='tilt',
                      default=0.0,
                      help='Desired tilt angle (deg)')
    parser.add_option('--rot',
                      dest='rot',
                      default=90.0,
                      help='Desired rotation angle (deg)')
    parser.add_option('--drift',
                      dest='drift',
                      default=0.0,
                      help='Desired drift angle (deg)')
    parser.add_option('--hdg',
                      dest='hdg',
                      default=45.0,
                      help='Desired hdg angle (deg)')
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

    global _el, _az
    global _tilt, _rot, _hdg, _drift
    global _pitchMin, _pitchMax, _pitchDelta
    global _rollMin, _rollMax, _rollDelta

    _el = float(options.el)
    _az = float(options.az)

    _tilt = float(options.tilt)
    _rot = float(options.rot)
    _hdg = float(options.hdg)
    _drift = float(options.drift)

    _pitchMin = float(options.pitchMin)
    _pitchMax = float(options.pitchMax)
    _pitchDelta = float(options.pitchDelta)

    _rollMin = float(options.rollMin)
    _rollMax = float(options.rollMax)
    _rollDelta = float(options.rollDelta)

    print("###############################", file=sys.stderr)
    print("#  el: ", _el, file=sys.stderr)
    print("#  az: ", _az, file=sys.stderr)
    print("#  tilt: ", _tilt, file=sys.stderr)
    print("#  rot: ", _rot, file=sys.stderr)
    print("#  drift: ", _drift, file=sys.stderr)
    print("#  hdg: ", _hdg, file=sys.stderr)
    print("#  pitchMin: ", _pitchMin, file=sys.stderr)
    print("#  pitchMax: ", _pitchMax, file=sys.stderr)
    print("#  pitchDelta: ", _pitchDelta, file=sys.stderr)
    print("#  rollMin: ", _rollMin, file=sys.stderr)
    print("#  rollMax: ", _rollMax, file=sys.stderr)
    print("#  rollDelta: ", _rollDelta, file=sys.stderr)

    # compute angles, and write them out

    # computeAngles()
    computeAngles2()

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

    print("#############################################", file=sys.stderr)
    print("# ",
          '{:>10} '.format('pitch'),
          '{:>10} '.format('roll'),
          '{:>10} '.format('el'),
          '{:>10} '.format('az'),
          '{:>10} '.format('rot'),
          '{:>10} '.format('tilt'),
          '{:>10} '.format('reflRot'),
          '{:>10} '.format('reflTilt'),
          file=sys.stderr)

    for pitch in np.arange(_pitchMin, _pitchMax, _pitchDelta):

        for roll in np.arange(_rollMin, _rollMax, _rollDelta):

            el, az = computeAzElYPrime(pitch, roll, _hdg, _rot, _tilt)

            rot, tilt = computeRotTiltYPrime(pitch, roll, _hdg, el, az)

            reflRot, reflTilt = computeReflAngles(pitch, roll, _drift, _rot, _tilt)

            print("  ",
                  '{:10.4f} '.format(pitch),
                  '{:10.4f} '.format(roll),
                  '{:10.4f} '.format(el),
                  '{:10.4f} '.format(az),
                  '{:10.4f} '.format(rot),
                  '{:10.4f} '.format(tilt),
                  '{:10.4f} '.format(reflRot),
                  '{:10.4f} '.format(reflTilt),
                  file=sys.stderr)

########################################################################
# Compute the angles, and write them out

def computeAngles2():

    print("#############################################", file=sys.stderr)
    print("# ",
          '{:>10} '.format('pitch'),
          '{:>10} '.format('roll'),
          '{:>10} '.format('rot'),
          '{:>10} '.format('tilt'),
          '{:>10} '.format('reflRot'),
          '{:>10} '.format('reflTilt'),
          file=sys.stderr)

    for pitch in np.arange(_pitchMin, _pitchMax, _pitchDelta):

        for roll in np.arange(_rollMin, _rollMax, _rollDelta):

            rot, tilt = computeRotTiltYPrime(pitch, roll, _hdg, _el, _az)

            reflRot, reflTilt = computeReflAngles(pitch, roll, _drift,
                                                  90.0 - _el, _az)

            print("  ",
                  '{:10.4f} '.format(pitch),
                  '{:10.4f} '.format(roll),
                  '{:10.4f} '.format(rot),
                  '{:10.4f} '.format(tilt),
                  '{:10.4f} '.format(reflRot),
                  '{:10.4f} '.format(reflTilt),
                  file=sys.stderr)

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
    len_a = math.sqrt(x_a * x_a + y_a * y_a + z_a * z_a)

    # compute matrix elements after multiplication
    # for 3 axis transformation

    mf11 = cosHdg * cosRoll + sinHdg * sinPitch * sinRoll
    mf12 = sinHdg * cosPitch
    mf13 = cosHdg * sinRoll - sinHdg * sinPitch * cosRoll

    mf21 = -sinHdg * cosRoll + cosHdg * sinPitch * sinRoll
    mf22 = cosHdg * cosPitch
    mf23 = -sinHdg * sinRoll - cosHdg * sinPitch * cosRoll

    mf31 = -cosPitch * sinRoll
    mf32 = sinPitch
    mf33 = cosPitch * cosRoll

    # Compute unit vector in earth coords

    xx = mf11 * x_a + mf12 * y_a + mf13 * z_a
    yy = mf21 * x_a + mf22 * y_a + mf23 * z_a
    zz = mf31 * x_a + mf32 * y_a + mf33 * z_a
    len = math.sqrt(xx * xx + yy * yy + zz * zz)

    # compute az and el

    az = math.degrees(math.atan2(xx, yy))
    el = math.degrees(math.asin(zz))

    if (options.debug):
        print("############# computeAzElYPrime ##############", file=sys.stderr)
        print("#  pitch: ", '{:10.4f} '.format(pitch), file=sys.stderr)
        print("#  roll : ", '{:10.4f} '.format(roll), file=sys.stderr)
        print("#  hdg  : ", '{:10.4f} '.format(hdg), file=sys.stderr)
        print("#  rot  : ", '{:10.4f} '.format(rot), file=sys.stderr)
        print("#  tilt : ", '{:10.4f} '.format(tilt), file=sys.stderr)

        print("#  x_a  : ", '{:10.4f} '.format(x_a), file=sys.stderr)
        print("#  y_a  : ", '{:10.4f} '.format(y_a), file=sys.stderr)
        print("#  z_a  : ", '{:10.4f} '.format(z_a), file=sys.stderr)
        print("#  len_a: ", '{:10.4f} '.format(len_a), file=sys.stderr)

        print("#  xx   : ", '{:10.4f} '.format(xx), file=sys.stderr)
        print("#  yy   : ", '{:10.4f} '.format(yy), file=sys.stderr)
        print("#  zz   : ", '{:10.4f} '.format(zz), file=sys.stderr)
        print("#  len  : ", '{:10.4f} '.format(len), file=sys.stderr)

        print("#  az   : ", '{:10.4f} '.format(az), file=sys.stderr)
        print("#  el   : ", '{:10.4f} '.format(el), file=sys.stderr)

    return (el, az)

########################################################################
# compute (elevation, azimuth) from rotation and tilt, attitude
# For a Y-prime radar e.g. HCR

def computeRotTiltYPrime(pitch, roll, hdg, el, az):

    # precompute sin/cos
    
    sinPitch = math.sin(math.radians(pitch))
    cosPitch = math.cos(math.radians(pitch))

    sinRoll = math.sin(math.radians(roll))
    cosRoll = math.cos(math.radians(roll))

    sinHdg = math.sin(math.radians(hdg))
    cosHdg = math.cos(math.radians(hdg))

    sinEl = math.sin(math.radians(el))
    cosEl = math.cos(math.radians(el))
    
    sinAz = math.sin(math.radians(az))
    cosAz = math.cos(math.radians(az))

    # compute unit vector relative to aircraft

    xx = sinAz * cosEl
    yy = cosAz * cosEl
    zz = sinEl
    len = math.sqrt(xx * xx + yy * yy + zz * zz)

    # compute matrix elements after multiplication
    # for 3 axis transformation

    mr11 = cosRoll * cosHdg + sinRoll * sinPitch * sinHdg
    mr12 = -cosRoll * sinHdg + sinRoll * sinPitch * cosHdg
    mr13 = -sinRoll * cosPitch

    mr21 = cosPitch * sinHdg
    mr22 = cosPitch * cosHdg
    mr23 = sinPitch

    mr31 = sinRoll * cosHdg - cosRoll * sinPitch * sinHdg
    mr32 = -sinRoll * sinHdg - cosRoll * sinPitch * cosHdg
    mr33 = cosRoll * cosPitch

    # Compute unit vector in earth coords

    x_a = mr11 * xx + mr12 * yy + mr13 * zz
    y_a = mr21 * xx + mr22 * yy + mr23 * zz
    z_a = mr31 * xx + mr32 * yy + mr33 * zz
    len_a = math.sqrt(x_a * x_a + y_a * y_a + z_a * z_a)

    # compute rot and tilt

    tilt = math.degrees(math.asin(y_a))
    rot = math.degrees(math.atan2(x_a, z_a))

    if (options.debug):
        print("############# computeRotTiltYPrime ##############", file=sys.stderr)
        print("#  pitch: ", '{:10.4f} '.format(pitch), file=sys.stderr)
        print("#  roll : ", '{:10.4f} '.format(roll), file=sys.stderr)
        print("#  hdg  : ", '{:10.4f} '.format(hdg), file=sys.stderr)
        print("#  rot  : ", '{:10.4f} '.format(rot), file=sys.stderr)
        print("#  tilt : ", '{:10.4f} '.format(tilt), file=sys.stderr)

        print("#  xx   : ", '{:10.4f} '.format(xx), file=sys.stderr)
        print("#  yy   : ", '{:10.4f} '.format(yy), file=sys.stderr)
        print("#  zz   : ", '{:10.4f} '.format(zz), file=sys.stderr)
        print("#  len  : ", '{:10.4f} '.format(len), file=sys.stderr)

        print("#  x_a  : ", '{:10.4f} '.format(x_a), file=sys.stderr)
        print("#  y_a  : ", '{:10.4f} '.format(y_a), file=sys.stderr)
        print("#  z_a  : ", '{:10.4f} '.format(z_a), file=sys.stderr)
        print("#  len_a: ", '{:10.4f} '.format(len_a), file=sys.stderr)

        print("#  az   : ", '{:10.4f} '.format(az), file=sys.stderr)
        print("#  el   : ", '{:10.4f} '.format(el), file=sys.stderr)

    return (rot, tilt)

########################################################################
# compute reflector angles
# pass in georef and desired rotation and tilt
# Returns rotation and tilt for the reflector
# This is from MotionControl.cpp.

def computeReflAngles(pitch, roll, drift, rot, tilt):

    sinPitch = math.sin(math.radians(pitch))
    cosPitch = math.cos(math.radians(pitch))

    sinRoll = math.sin(math.radians(roll))
    cosRoll = math.cos(math.radians(roll))

    sinDrift = math.sin(math.radians(drift))
    cosDrift = math.cos(math.radians(drift))

    # Track relative coordinates - desired beam position

    sinRot = math.sin(math.radians(rot))
    cosRot = math.cos(math.radians(rot))

    sinTilt = math.sin(math.radians(tilt))
    cosTilt = math.cos(math.radians(tilt))

    # Convert to track relative Cartesian coordinates

    xx = cosTilt * sinRot
    yy = sinTilt
    zz = cosTilt * cosRot
    len = math.sqrt(xx * xx + yy * yy + zz * zz)

    # Convert to pod relative Cartesian coordinates - adjusted beam position

    mr11 = cosDrift * cosRoll + sinDrift * sinPitch * sinRoll
    mr12 = -sinDrift * cosRoll + cosDrift * sinPitch * sinRoll
    mr13 = -cosPitch * sinRoll

    mr21 = sinDrift * cosPitch
    mr22 = cosDrift * cosPitch
    mr23 = sinPitch

    mr31 = cosDrift * sinRoll - sinDrift * sinPitch * cosRoll
    mr32 = -sinDrift * sinRoll - cosDrift * sinPitch * cosRoll
    mr33 = cosPitch * cosRoll

    x_a = mr11 * xx + mr12 * yy + mr13 * zz
    y_a = mr21 * xx + mr22 * yy + mr23 * zz
    z_a = mr31 * xx + mr32 * yy + mr33 * zz
    len_a = math.sqrt(x_a * x_a + y_a * y_a + z_a * z_a)

    # x_a = \
    #       xx * (cosDrift * cosRoll - sinDrift * sinPitch * sinRoll) + \
    #       yy * (sinDrift + cosDrift * sinPitch * sinRoll) + \
    #       -zz * cosPitch * sinRoll

    # y_a = \
    #       -xx * sinDrift * cosPitch + \
    #       yy * cosDrift * cosPitch + \
    #       zz * sinPitch

    # z_a = \
    #       xx * (cosDrift * sinRoll + sinDrift * sinPitch * sinRoll) + \
    #       yy * (sinDrift * sinRoll - cosDrift * sinPitch * cosRoll) + \
    #       zz * cosPitch * cosRoll

    # len_a = math.sqrt(x_a * x_a + y_a * y_a + z_a * z_a)

    # Convert from pod relative Cartesian coordinates to polar coordinates
    # and save the adjusted rotation and tilt angles.
    # devide by 2 because it is a reflector

    reflTilt = math.degrees(math.asin(y_a)) / 2.0

    # KLUGE: The algorithm above isn't really right. As compensation for now,
    # just change the sign of the corrected tilt if the desired rotation angle
    # rots downward, i.e., if its cosine is less than zero.

    if (cosRot < 0.0):
        reflTilt = reflTilt * -1

    # rotation needs to be in the range of 0-360

    reflRot = math.degrees(math.atan2(x_a, z_a))

    #if (reflRot < 0):
    #    reflRot = reflRot + 360.0

    if (options.debug):
        print("############# computeRefAngles ##############", file=sys.stderr)
        print("#  pitch: ", '{:10.4f} '.format(pitch), file=sys.stderr)
        print("#  roll : ", '{:10.4f} '.format(roll), file=sys.stderr)
        print("#  drift: ", '{:10.4f} '.format(drift), file=sys.stderr)
        print("#  rot  : ", '{:10.4f} '.format(rot), file=sys.stderr)
        print("#  tilt : ", '{:10.4f} '.format(tilt), file=sys.stderr)

        print("#  xx   : ", '{:10.4f} '.format(xx), file=sys.stderr)
        print("#  yy   : ", '{:10.4f} '.format(yy), file=sys.stderr)
        print("#  zz   : ", '{:10.4f} '.format(zz), file=sys.stderr)
        print("#  len  : ", '{:10.4f} '.format(len), file=sys.stderr)

        print("#  x_a  : ", '{:10.4f} '.format(x_a), file=sys.stderr)
        print("#  y_a  : ", '{:10.4f} '.format(y_a), file=sys.stderr)
        print("#  z_a  : ", '{:10.4f} '.format(z_a), file=sys.stderr)
        print("#  len_a: ", '{:10.4f} '.format(len_a), file=sys.stderr)

        print("#  reflRot  : ", '{:10.4f} '.format(reflRot), file=sys.stderr)
        print("#  reflTilt : ", '{:10.4f} '.format(reflTilt), file=sys.stderr)

    return (reflRot, reflTilt)

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

