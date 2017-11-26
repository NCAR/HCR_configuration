#!/bin/sh
#
# Send 20 IWG1 packets to multicast address 239.0.0.1 on port 7071.
# This packet has lat/lon/alt 39.9/-105.1/1700, heading 89.9, and tas 0.0,
# which will make the Cmigits class perform a stationary initialization.
#
NPKTS=20
MCAST_GROUP="239.0.0.10"
MCAST_PORT=7071
MCAST_TTL=2	# set this bigger if the multicast needs to make more hops

echo -n "Sending $nPkts IWG1 packets"
for ((i = 0 ; i < $NPKTS ; i++ )); do
    echo -n "."
    echo "IWG1,20131003T145500,39.9,-105.1,1700.,,4259.16,4238.06,0.0,0.0,209.665,0.163316,-0.0257365,89.9,189.041,4.14272,2.9284,1.11889,0.413966,2.94089,5.6283,-9.3279,14.2788,614.518,73.1688,856.621,8.34496,83.4357,0.10949,0.81966,0.751137,," | socat STDIO UDP-DATAGRAM:$MCAST_GROUP:$MCAST_PORT,ip-multicast-ttl=$MCAST_TTL
    sleep 1
done
echo "done"
