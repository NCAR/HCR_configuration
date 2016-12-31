## Installing the HCR configuration

### Checking out the HCR repos

Make sure you have checked out the HCR run-time configuration.

See [README_HCR_CHECKOUT.md](./README_HCR_CHECKOUT.md)

### Make the top-level data directory

Create the following directory, owned by the hcr user:

```
  /data/hcr
```

### Run the configuration script

Go to the system scripts directory:

```
  cd ~/git/HCR_configuration/projDir/system/scripts
```

Run the script:

```
  ./configureHost.py
```

There are 2 host types:

```
  archiver
  drx
```

You need to specify the relevant host type.

### Check the setup and links:

In the home directory, you should see the following links:

```
  .Xdefaults -> /home/hcr/git/HCR_configuration/projDir/system/dotfiles/Xdefaults
  .cshrc -> /home/hcr/git/HCR_configuration/projDir/system/dotfiles/cshrc
  .emacs -> /home/hcr/git/HCR_configuration/projDir/system/dotfiles/emacs
  projDir -> /home/hcr/git/HCR_configuration/projDir
```

On the **archiver** host you should see:


```

  cd ~/projDir
  data -> /data/hcr/data.archiver
  logs -> data/logs

  cd ~/projDir/control
  data_list -> data_list.archiver
  proc_list -> proc_list.archiver

```

On the **drx** host you should see:


```

  cd ~/projDir
  data -> /data/hcr/data.drx
  logs -> data/logs

  cd ~/projDir/control
  data_list -> data_list.drx
  proc_list -> proc_list.drx

```

### Running the system

To start the system:

```
  source ~/.cshrc
  start_all
```

To monitor the system

```

  # see running processes
  ppm

  # check all processes are running
  pcheck

  # check data flow
  pdm

```

To stop the system:

```
  stop_all
```


