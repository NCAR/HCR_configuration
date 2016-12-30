#!/bin/sh
mcsFile=~hcr/workspace/sd3c/firmware/bitstream/sd3c-ddc8.mcs
echo "Loading ${mcsFile##*/} to the Pentek"
~hcr/HCR_instrument/src/ReadyFlow7142_428/linux/1.0/x86_64/examples/fpgaload.out $mcsFile
