#!/bin/sh

# Copy all data except time-series from external HCR* disks to the NAS
# synced_data/hcr directory. Data here get synced back to Boulder 
# automatically.
rsync -av --exclude time_series --exclude today --exclude yesterday \
    /run/media/hcr/HCR*/hcr/archive/CSET/* \
    /nas/data/synced_data/hcr/

# Copy all data except time-series from external HCR* disks to the NAS
# scr_data/hcr directory. This directory is used by local data displays.
rsync -av --exclude time_series --exclude today --exclude yesterday \
    /run/media/hcr/HCR*/hcr/archive/CSET/* \
    /nas/data/scr_data/hcr/
