#!/usr/bin/env python

#===========================================================================
#
# Adjust HCR velocity for antenna angles
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
                      default=10.0,
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

    global _el, _az
    global _tilt, _rot, _hdg, _drift
    global _pitchMin, _pitchMax, _pitchDelta
    global _rollMin, _rollMax, _rollDelta

    _el = float(options.el)
    _az = float(options.az)

    _tilt = float(options.tilt)
    _rot = float(options.rot)
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
    print("#  pitchMin: ", _pitchMin, file=sys.stderr)
    print("#  pitchMax: ", _pitchMax, file=sys.stderr)
    print("#  pitchDelta: ", _pitchDelta, file=sys.stderr)
    print("#  rollMin: ", _rollMin, file=sys.stderr)
    print("#  rollMax: ", _rollMax, file=sys.stderr)
    print("#  rollDelta: ", _rollDelta, file=sys.stderr)

    # compute corrections for series of angles

    computeCorrections()

    sys.exit(0)
    
########################################################################
# Compute the angles, and write them out

def computeCorrections():
    
    print("#############################################", file=sys.stderr)
    print("# ",
          '{:>10} '.format('pitch'),
          '{:>10} '.format('roll'),
          '{:>10} '.format('drift'),
          '{:>10} '.format('rot'),
          '{:>10} '.format('tilt'),
          '{:>10} '.format('velCorr'),
          file=sys.stderr)

    for pitch in np.arange(_pitchMin, _pitchMax, _pitchDelta):

        for roll in np.arange(_rollMin, _rollMax, _rollDelta):
            
            velCorr = computeVelCorr(pitch, roll, _drift, _rot, _tilt,
                                     100.0, 50.0, 5.0, 5.0)

            print("  ",
                  '{:10.4f} '.format(pitch),
                  '{:10.4f} '.format(roll),
                  '{:10.4f} '.format(_drift),
                  '{:10.4f} '.format(_rot),
                  '{:10.4f} '.format(_tilt),
                  '{:10.4f} '.format(velCorr),
                  file=sys.stderr)

########################################################################
# compute vel correction from rotation and tilt, and attitude

def computeVelCorr(pitch, roll, drift, rot, tilt,
                   ns_vel, ew_vel, vert_vel,
                   vel_meas):

    # precompute sin/cos
    
    sinPitch = math.sin(math.radians(pitch))
    cosPitch = math.cos(math.radians(pitch))

    sinRoll = math.sin(math.radians(roll))
    cosRoll = math.cos(math.radians(roll))

    sinDrift = math.sin(math.radians(drift))
    cosDrift = math.cos(math.radians(drift))

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
    
    mf11 = cosDrift * cosRoll + -sinDrift * sinPitch * sinRoll
    mf12 = -sinDrift * cosPitch
    mf13 = cosDrift * sinRoll + sinDrift * sinPitch * cosRoll
    
    mf21 = sinDrift * cosRoll + cosDrift * sinPitch * sinRoll
    mf22 = cosDrift * cosPitch
    mf23 = sinDrift * sinRoll - cosDrift * sinPitch * cosRoll

    mf31 = -cosPitch * sinRoll
    mf32 = sinPitch
    mf33 = cosPitch * cosRoll

    # Compute unit vector in track-centered coords
    # Platform is horizontal, and yy aligns with the track

    xx = mf11 * x_a + mf12 * y_a + mf13 * z_a
    yy = mf21 * x_a + mf22 * y_a + mf23 * z_a
    zz = mf31 * x_a + mf32 * y_a + mf33 * z_a
    len_t = math.sqrt(xx * xx + yy * yy + zz * zz)

    #
    # Since this is along-track we can say:
    #
    #   xx: there is no platform velocity in xx dirn
    #   yy: all horizontal motion is in yy dirn
    #   zz: vertical motion is in zz dirn

    # ground vel is towards platform
    # so for positive yy, the platform motion will induce a negative bias
    # so the correction must be positive for positive yy

    ground_speed = math.sqrt(ns_vel * ns_vel + ew_vel * ew_vel)
    horiz_corr = yy * ground_speed

    # for positive vert_vel, the measured velocity when pointing up
    # will be negative, since the particles appear to approach the plaform
    # so the correction is positive for positive vert_vel

    vert_corr = zz * vert_vel
    
    # since the two corrections are othogonal, they can be summed

    vel_corr = vel_meas + horiz_corr + vert_corr

    if (options.debug):
        print("############# computeAzElYPrime ##############", file=sys.stderr)
        print("#  pitch   : ", '{:10.4f} '.format(pitch), file=sys.stderr)
        print("#  roll    : ", '{:10.4f} '.format(roll), file=sys.stderr)
        print("#  drift   : ", '{:10.4f} '.format(drift), file=sys.stderr)
        print("#  rot     : ", '{:10.4f} '.format(rot), file=sys.stderr)
        print("#  tilt    : ", '{:10.4f} '.format(tilt), file=sys.stderr)

        print("#  ns_vel  : ", '{:10.4f} '.format(ns_vel), file=sys.stderr)
        print("#  ew_vel  : ", '{:10.4f} '.format(ew_vel), file=sys.stderr)
        print("#  vert_vel: ", '{:10.4f} '.format(vert_vel), file=sys.stderr)

        print("#  vel_meas: ", '{:10.4f} '.format(vel_meas), file=sys.stderr)

        print("#  x_a     : ", '{:10.4f} '.format(x_a), file=sys.stderr)
        print("#  y_a     : ", '{:10.4f} '.format(y_a), file=sys.stderr)
        print("#  z_a     : ", '{:10.4f} '.format(z_a), file=sys.stderr)
        print("#  len_a   : ", '{:10.4f} '.format(len_a), file=sys.stderr)

        print("#  xx      : ", '{:10.4f} '.format(xx), file=sys.stderr)
        print("#  yy      : ", '{:10.4f} '.format(yy), file=sys.stderr)
        print("#  zz      : ", '{:10.4f} '.format(zz), file=sys.stderr)
        print("#  len_t   : ", '{:10.4f} '.format(len_t), file=sys.stderr)

        print("#  ground_speed : ", '{:10.4f} '.format(ground_speed), file=sys.stderr)
        print("#  horiz_corr   : ", '{:10.4f} '.format(horiz_corr), file=sys.stderr)
        print("#  vert_corr    : ", '{:10.4f} '.format(vert_corr), file=sys.stderr)
        
        print("#  vel_corr     : ", '{:10.4f} '.format(vel_corr), file=sys.stderr)

    return vel_corr


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

