# Running the HCR system

## To start processes on the rds (drx) host


```
  start_all
  start_hcrdrx.ops
```

To monitor:

```

  # see running processes
  ppm

  # check all processes are running
  pcheck

  # check data flow
  pdm

```

## To start the system on the archiver host

```
  start_all
```

To monitor:

```

  # see running processes
  ppm

  # check all processes are running
  pcheck

  # check data flow
  pdm

```

NOTE: TerrainHtServer takes a while to start up and register with procmap.

## Running in the lab

Turn on the power switches for the pod and the nosecone.

After the system has come up, run the following on the archiver:

```
  send_iwg1_packet.sh
```

## To stop the system

On each host run:

```
  stop_all
```


