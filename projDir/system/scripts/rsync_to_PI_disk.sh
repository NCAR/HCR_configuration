#!/bin/bash
#
# This script rsyncs HCR data from an HCR disk to a PI-provided disk for OTREC
# To sync data from HCR source disk labeled RSF0012, the command is just:
#
#    rsync_to_PI_disk.sh RSF0012
#
# The script finds the destination PI disk automatically. The PI disks for
# OTREC don't have filesystem labels, so they are mounted with names based
# on their UUIDs, e.g.,
#
#     /run/media/hcr/65f8760d-b8d9-435a-bdf4-2f76e13380c8
#
# We look for a mounted disk with a name of this form, and if we find exactly
# one, we use it as the destination disk. Otherwise we print an error and exit.

if [ $# -ne 1 ] || [ "$1" == "-h" ]; then
    echo "This script rsyncs OTREC data from an HCR disk to a PI-provided"
    echo "disk. The script will try to find the PI disk name automatically."
    echo "Example to rsync data from source HCR disk RSF0012:"
    echo " "
    echo "    $ $0 RSF0012"
    echo " "
    echo "Usage: $0 <HCR_source_disk>"
    exit 1
fi

# Top mount directory for the user's USB disks
usb_disk_prefix="/run/media/${USER}"

# Source disk
sourceDisk="${usb_disk_prefix}/$1"

# Regex to match a 36-character UUID line containing only 0-9, a-f, or '-'
uuid_regex='^[0-9a-f-]{36}$'

# Find any USB disks mounted by UUID
usb_uuids=`ls ${usb_disk_prefix} | grep -E $uuid_regex`
match_count=`echo $usb_uuids | wc -w`
if [ $match_count == 0 ]; then
    echo ""
    echo "Exiting because no destination PI disk (named by UUID) was found!"
    exit 1
elif [ $match_count != 1 ]; then
    echo "Exiting because $match_count PI disks (named by UUID) were found. Unmount all but one!"
    exit 1
fi

# PI disk path
piDisk=${usb_disk_prefix}/${usb_uuids}

# If necessary, make a top-level "rsf" directory on the PI disk, with world
# write permission
#if [ ! -x ${piDisk}/rsf || ! -w ${piDisk}/rsf ]; then
if [ ! -w ${piDisk}/rsf ]; then
    echo "Adding 'rsf' dir with write permission"
    echo "    to PI disk ${piDisk}"
    echo "ENTER YOUR PASSWORD BELOW FOR sudo PERMISSION"
    sudo -k   # force sudo timeout so the following always requests a password
    sudo mkdir -p -m 777 ${piDisk}/rsf || exit 1
fi

# Clone the HCR source disk to the destination PI disk using rsync
# The -O flag avoids an error message about setting time on the top level
# destination directory.
echo "Cloning ${sourceDisk} -> ${piDisk}"
rsync -O -a --info progress2 ${sourceDisk}/rsf/ ${piDisk}/rsf/

exit 0
