#!/bin/sh -x
mcsFile=~hcr/workspace/sd3c/firmware/bitstream/sd3c-ddc8.mcs

echo "Loading ${mcsFile##*/} to the Pentek"
ReadyFlowDir=~hcr/git/HCR_instrument/src/ReadyFlow7142_428_WinDriver1150/linux/1.0/x86_64
LD_LIBRARY_PATH=$ReadyFlowDir/lib
$ReadyFlowDir/examples/fpgaload.out $mcsFile
