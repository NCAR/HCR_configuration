# Configuring the archiver host

The software on both the instrument and archiver computers run under Centos LINUX.

The user is **hcr**.
Set the hsrl user login shell to **tcsh**.
Make sure **python** is installed.

To install the LROSE software, see:

  [LROSE core](https://ncar.github.io/lrose-core)

To configure the runtime environment on the HSRL archiver host, first check out hsrl_configuration:

```bash
  mkdir -p ~/git
  cd ~/git
  git clone https://github.com/ncar/HCR_configuration
```

Make the data directory, and give it the correct permissions:

```bash
  sudo mkdir -p /data/hcr
  sudo chown -R hcr /data/hcr
```

Go to the systems/scripts directory, and run the configuration script:

```bash
  cd ~/git/HCR_configuration/projDir/system/scripts
  ./configureHost.py
```

Source the .cshrc file

```csh
  source ~/.cshrc
```

Go to the projects directory:

```csh
  cd ~/projDir
```
