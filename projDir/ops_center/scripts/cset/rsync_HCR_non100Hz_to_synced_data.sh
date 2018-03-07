#!/bin/sh

# Copy 10 Hz CFRadial, logs, and spdb files from external HCR* disks to the NAS
# synced_data/hcr directory
rsync -av --exclude time_series --exclude 100hz \
    --exclude today --exclude yesterday \
    /run/media/hcr/HCR*/hcr/archive/CSET/* \
    /nas/data/synced_data/hcr/
