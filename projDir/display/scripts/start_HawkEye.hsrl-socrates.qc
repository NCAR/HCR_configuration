#! /bin/csh

# start BSCAN in archive mode

cd $PROJ_DIR/display/params

if ($#argv < 6) then
  HawkEye -params HawkEye.hsrl-socrates.qc -start_time "2018 02 04 02 30 00" >& /dev/null &
else
  HawkEye -params HawkEye.hsrl-socrates.qc -start_time "$1 $2 $3 $4 $5 $6" >& /dev/null &
endif

