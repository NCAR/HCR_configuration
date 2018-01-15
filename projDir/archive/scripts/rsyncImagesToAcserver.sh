#!/bin/sh
#
# Copy catalog image files to acserver to make them available to other machines
# on the plane.

# Get dates for today and yesterday in yyyymmdd form
todayDate=`date +%Y%m%d`
yesterdayDate=`date +%Y%m%d -d 'now - 1 day'`

# Rsync images from yesterday and today to acserver
rsync -av $PROJ_DIR/data/images/catalog/$yesterdayDate \
          ads@acserver:/var/www/html/hcr/images/
rsync -av $PROJ_DIR/data/images/catalog/$todayDate \
          ads@acserver:/var/www/html/hcr/images/

# Rsync HSRL images from yesterday and today to acserver
rsync -av $PROJ_DIR/data/hsrl/images/catalog/$yesterdayDate \
          ads@acserver:/var/www/html/hsrl/images/
rsync -av $PROJ_DIR/data/hsrl/images/catalog/$todayDate \
          ads@acserver:/var/www/html/hsrl/images/
