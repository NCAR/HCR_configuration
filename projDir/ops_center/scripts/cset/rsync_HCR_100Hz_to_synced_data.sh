#!/bin/sh

# Copy 100 Hz CFRadial files from external HCR* disks to the NAS
# synced_data/hcr directory
rsync -av --exclude 10hz \
    /run/media/hcr/HCR*/hcr/archive/CSET/cfradial \
    /nas/data/synced_data/hcr/
