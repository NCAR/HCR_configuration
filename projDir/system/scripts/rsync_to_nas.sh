#!/bin/bash
#
# This script rsyncs HCR and HSRL data from an HCR disk to the
# OTREC ops center server. To sync data from disk HCR39, the
# command is just:
#
#    rsync_to_nas.sh RSF0012
#
if [ $# -ne 1 ] || [ "$1" == "-h" ]; then
    echo "This script rsyncs data from an HCR disk to the OTREC"
    echo "ops center server. E.g., to rsync data from disk RSF0012:"
    echo " "
    echo "    $ $0 RSF0012"
    echo " "
    echo "Usage: $0 <disk_name>"
    exit 1
fi

hcrDisk=$1

# rsync HCR data which will automatically be synced to Boulder
rsync -av --info=flist1,stats0 --exclude hsrl --exclude time_series \
    --exclude 100hz --exclude yesterday --exclude today --exclude '*.tmp' \
    /run/media/hcr/$hcrDisk/rsf/archive/otrec/* \
    /nas/data/data/OTREC/hcr_synced/

if [ ! $? -eq 0 ]; then
    echo "HCR RSYNC FAILED!"
    exit 1
fi

# rsync HCR 100 Hz data (NOT automatically synced to Boulder for
# bandwidth reasons)
rsync -av --info=flist1,stats0 --exclude 10hz --exclude '*.tmp' \
    /run/media/hcr/$hcrDisk/rsf/archive/otrec/cfradial \
    /nas/data/data/OTREC/hcr_100hz/

if [ ! $? -eq 0 ]; then
    echo "HCR 100hz RSYNC FAILED!"
    exit 1
fi

## then rsync HSRL data which will automatically be synced to Boulder
#rsync -av --info=flist1,stats0 \
#    /run/media/hcr/$hcrDisk/rsf/archive/otrec/hsrl/* \
#    /nas/data/data/OTREC/hsrl_synced/
#
#if [ ! $? -eq 0 ]; then
#    echo "HSRL RSYNC FAILED!"
#    exit 1
#fi

#echo "HCR and HSRL data syncs were successful"
echo "HCR data syncs were successful"
exit 0


