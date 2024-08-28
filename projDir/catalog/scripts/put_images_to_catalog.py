#!/usr/bin/env python

#===========================================================================
#
# Put image data to catalog
#
#===========================================================================



import os
import sys
from stat import *
import time
import datetime
from datetime import timedelta
import string
import logging
import paramiko
import subprocess
from optparse import OptionParser

# Set up our basic logging configuration, writing to stdout with a simple format
logging.basicConfig(
    level=logging.INFO,
    format='%(levelname)s - %(message)s',
    stream=sys.stdout
)

# Create our global logger
logger = logging.getLogger()

def main():

    global options

    # parse the command line

    parseArgs()

    # Set log level to DEBUG if requested
    if (options.debug):
        print("ENABLING DEBUG")
        logger.setLevel(logging.DEBUG)

    # Log the parsed options
    logger.debug("Options:")
    logger.debug("  debug/verbose? " + str(options.debug))
    logger.debug("  validTime: " + str(options.validTime))
    logger.debug("  imageDir: " + options.imageDir)
    logger.debug("  relDataPath: " + options.relDataPath)
    logger.debug("  fileName: " + options.fileName)
    logger.debug("  ftpServer: " + options.ftpServer)
    logger.debug("  targetDir: " + options.targetDir)
    logger.debug("  category: " + options.category)
    logger.debug("  platform: " + options.platform)
    logger.debug("  removeAfterCopy: " + str(options.removeAfterCopy))

    # load our SFTP user name and password into global variables sftpUser and
    # sftpPassword
    loadSFTPUserSecrets()

    # initialize

    logger.debug("==================================================================")
    logger.debug("BEGIN: put_images_to_catalog.py at " + str(datetime.datetime.now()))
    logger.debug("==================================================================")

    #   compute valid time string

    validTime = time.gmtime(int(options.validTime))
    year = int(validTime[0])
    month = int(validTime[1])
    day = int(validTime[2])
    hour = int(validTime[3])
    minute = int(validTime[4])
    sec = int(validTime[5])
    validDayStr = "%.4d%.2d%.2d" % (year, month, day)
    validTimeStr = "%.4d%.2d%.2d%.2d%.2d" % (year, month, day, hour, minute)
    dateTimeStr = "%.4d%.2d%.2d-%.2d%.2d%.2d" % (year, month, day, hour, minute, sec)

    # compute full path of image

    fullFilePath = options.imageDir
    fullFilePath = os.path.join(fullFilePath, options.relDataPath);

    # extract the platform and product from the file name.
    # The image files are named like:
    #   <category>.<platform>.<time>.<field>.png.

    file_tokens = options.fileName.split(".")
    logger.debug("filename toks: ")
    logger.debug(file_tokens)

    if len(file_tokens) != 5:
        logger.error("Invalid fileName: %s" % ( options.fileName))
        sys.exit(1)

    # category

    category = file_tokens[0]
    if (options.category.find("NONE") < 0):
        category = options.category

    # platform

    platform = file_tokens[1]
    if (options.platform.find("NONE") < 0):
        platform = options.platform

    # platform

    timeStr = file_tokens[2]

    # field

    field_name = file_tokens[3]

    # time - use from latest data file
    # time_tokens = file_tokens[4].split(".")
    # validTimeStr = time_tokens[0];

    # compute catalog file name

    catalogName = (category + "." + platform + "." +
                   validTimeStr + "." +
                   field_name + "." + "png")

    logger.debug("catalogName: %s" % catalogName)

    # put the image file

    putFile(fullFilePath, catalogName)

    # optionally remove the file

    if (options.removeAfterCopy == True):
        logger.debug("Removing file: ", fullFilePath)
        os.remove(fullFilePath)

    # let the user know we are done

    logger.debug("==================================================================")
    logger.debug("END: put_images_to_catalog.py at " + str(datetime.datetime.now()))
    logger.debug("==================================================================")

    sys.exit(0)

########################################################################
# Load SFTP username and password from our SFTP user secrets file into
# global variables sftpUser and sftpPassword.
#
# The SFTP username and password are now kept in read-protected file
# /home/hcr/.ssh/sftp_user_secrets.py. Its contents are simple:
#
#       sftpUser = "<username>"
#       sftpPassword = "<password>"
#
# Just replace <username> and <password> with the desired values.
#
# Note that the file is treated as a Python source file, so it can
# contain comments.
#
# We use this simple secrets file so that the SFTP username and password do
# not get committed into the repository as part of this script.

def loadSFTPUserSecrets():

    sftp_secrets_file = "/home/hcr/.ssh/sftp_user_secrets.py"

    # Open our sftp user secrets file and execute its contents. On success,
    # global variables sftpUser and sftpPassword will hold the SFTP username
    # and password given in the file.
    global sftpUser; sftpUser = None
    global sftpPassword; sftpPassword=None

    try:
        with open(sftp_secrets_file) as secrets:
            exec(secrets.read(), globals())
    except Exception as e:
        logger.error("Failed to open or execute %s: %s" % (sftp_secrets_file, repr(e)))
        sys.exit(1)


########################################################################
# Put the specified file

def putFile(filePath, catalogName):

    logger.debug("Handling file: " + filePath)
    logger.debug("  catalogName: " + catalogName)

    # create tmp dir if necessary

    dataDir = os.environ['DATA_DIR']
    tmpDir = os.path.join(dataDir, "images/tmp")
    logger.debug("  tmpDir: " + tmpDir)
    if not os.path.exists(tmpDir):
        os.makedirs(tmpDir, 0o775)

    # copy the file to the tmp directory

    tmpPath = os.path.join(tmpDir, catalogName)
    cmd = "cp " + filePath + " " + tmpPath
    runCommand(cmd)

    # send the file to the catalog

    ftpFile(catalogName, tmpPath)

    # remove the tmp file

    os.remove(tmpPath)

    return 0

########################################################################
# Ftp the file

def ftpFile(fileName, filePath):

    targetDir = options.targetDir
    ftpServer = options.ftpServer

    # Handle SSH authentication with SFTP server

    sshClient = paramiko.SSHClient()
    sshClient.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    sshClient.connect(hostname=ftpServer, username=sftpUser, password=sftpPassword)

    # Now create a paramiko.SFTPClient instance using the open SSH connection

    sftpClient = sshClient.open_sftp()

    # if targetDir doesn't exist (i.e., if we can't list its contents), try to create it

    try:
        sftpClient.listdir(targetDir)
    except Exception:
        try:
            sftpClient.mkdir(targetDir)
            logger.info("Created remote target directory %s" % (targetDir))
        except Exception as e:
            logger.error("failed to create remote target directory %s: %s"
                         % (targetDir, e))
            sys.exit(1)

    # go to target dir

    logger.debug("sftp chdir to: " + targetDir)
    sftpClient.chdir(targetDir)

    # put the file

    logger.debug("putting file: " + filePath)

    sftpClient.put(filePath, fileName)

    # close ftp connection

    sftpClient.close()

    return

########################################################################
# Parse the command line

def parseArgs():

    global options

    # parse the command line

    usage = "usage: %prog [options]"
    parser = OptionParser(usage)

    # these options come from the ldata info file

    parser.add_option('--debug',
                      dest='debug', default=False,
                      action="store_true",
                      help='Set debugging on')
    parser.add_option('--verbose',
                      dest='debug', default=False,
                      action="store_true",
                      help='Synonym for --debug')
    parser.add_option('--unix_time',
                      dest='validTime',
                      default=0,
                      help='Valid time for image')
    parser.add_option('--full_path',
                      dest='imageDir',
                      default='unknown',
                      help='Full path of image file')
    parser.add_option('--file_name',
                      dest='fileName',
                      default='unknown',
                      help='Name of image file')
    parser.add_option('--rel_data_path',
                      dest='relDataPath',
                      default='unknown',
                      help='Relative path of image file')

    # these options are specific to the image type

    parser.add_option('--ftp_server',
                      dest='ftpServer',
                      default='catalog.eol.ucar.edu',
                      help='Target FTP server')
    parser.add_option('--target_dir',
                      dest='targetDir',
                      default='pub/incoming/catalog/winter',
                      help='Target directory on the FTP server')
    parser.add_option('--category',
                      dest='category',
                      default='NONE',
                      help='Category portion of the catalog file name')
    parser.add_option('--platform',
                      dest='platform',
                      default='NONE',
                      help='Platform portion of the catalog file name.  Overrides platform in image file name if specified.')
    parser.add_option('--remove_after_copy',
                      dest='removeAfterCopy', default=False,
                      action="store_true",
                      help='Remove files after copy to ftp server')

    (options, args) = parser.parse_args()

########################################################################
# Run a command in a shell, wait for it to complete

def runCommand(cmd):

    logger.debug("running local cmd:" + cmd)

    try:
        retcode = subprocess.call(cmd, shell=True)
        if retcode < 0:
            logger.error("Local shell command '%s' was terminated by signal %d"
                         % (cmd, -retcode))
        else:
            logger.debug("Local shell command '%s' returned code %d"
                         % (cmd, retcode))
    except OSError as e:
        logger.error("Executing local command '%s' failed: %s" % (cmd, e) )

########################################################################
# kick off main method

if __name__ == "__main__":

   main()
