#!/bin/sh
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/opt/aeros-5122M/deploy/lib
exec /opt/aeros-5122M/deploy/bin/aeros >& /dev/null &
