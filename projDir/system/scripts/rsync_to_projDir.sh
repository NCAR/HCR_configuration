#!/bin/bash
#
# This script rsyncs HCR and HSRL data from an HCR disk to the
# projDir on this host. To sync data from disk HCR39, the
# command is just:
#
#    rsync_to_projDir.sh HCR39
#
if [ $# -ne 1 ] || [ "$1" == "-h" ]; then
    echo "This script rsyncs data from an HCR disk to the projDir"
    echo "on the local machine. E.g., to rsync data from disk HCR39:"
    echo " "
    echo "    $ $0 HCR39"
    echo " "
    echo "Usage: $0 <disk_name>"
    exit 1
fi

hcrDisk=$1

# rsync data directories

rsync -av --info=flist1,stats0 --exclude time_series \
    --exclude 100hz --exclude yesterday --exclude today \
    /run/media/hcr/$hcrDisk/hcr/archive/socrates/* \
    ${HOME}/projDir/data

if [ ! $? -eq 0 ]; then
    echo "RSYNC FAILED!"
    exit 1
fi

echo "data sync was successful"
exit 0


