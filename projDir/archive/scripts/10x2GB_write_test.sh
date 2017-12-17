#!/bin/sh

# Destination directory
#destdir=/run/media/hcr/HCR54/foo
#destdir=/run/media/hcr/HCR105/foo
destdir=/run/media/mdtest/HCR000/foo
# Exit immediately on INT signal (i.e., ^C)
trap "exit" INT

# Source file
srcfile=/home/mdtest/Downloads/20150729_235134.868_999_141.iwrf_ts

# Source file size in bytes
srcsize=`du -sb $srcfile | cut -f1`

# Get srcfile into the memory cache for fast read access
echo -n "Getting the source file into memory cache: "
/bin/time -f "%e s" cat $srcfile > /dev/null

# Sync to disks before we start
sync

# Make 10 copies to the destination directory
copycount="10"
echo "Rsyncing $srcsize byte file into $destdir $copycount times..."
for num in `seq 1 $copycount`; do
  echo -n "Copy $num: "
  # Time rsync of the file to destdir. We sync afterwards to make sure the
  # copy to the dest disk is complete.
  /bin/time -f "%e s" sh -c "rsync $srcfile $destdir/test_$ && sync"
done
echo "Removing copied files from $destdir"
rm -f $destdir/test_*
