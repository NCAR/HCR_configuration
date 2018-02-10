#!/bin/sh

# Copy all data except time-series from external HCR* disks to the NAS
# scr_data/hcr directory
rsync -av --exclude time_series --exclude today --exclude yesterday \
    /run/media/hcr/HCR*/hcr/archive/CSET/* \
    /nas/data/scr_data/hcr/
