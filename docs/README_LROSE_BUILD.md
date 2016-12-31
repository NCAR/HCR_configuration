## Building using the LROSE development makefile system

### Ensure the required LINUX packages are installed.

See [README_LINUX_SETUP.md](./README_LINUX_SETUP.md)

### Checking out the HCR repos

See [README_HCR_CHECKOUT.md](./README_HCR_CHECKOUT.md)

### Build and install **netcdf support** in ~/lrose

```
  cd ~/git/lrose-netcdf
  ./build_and_install_netcdf -x ~/lrose
```

This will install netcdf in:

```
  ~/lrose/include
  ~/lrose/lib
  ~/lrose/bin
```

### Setting up the environment for the lrose build

The software development system at NCAR/RAL (formerly RAP) and NCAR/EOL makes use of a recursive makefile approach, using environment variables to identify the various directories used during the build.

If you have correctly installed the HCR configurations, and sourced the ~.,cshrc file
the build environment should already be set up.

If not, perform the actions in the following section.

##### Source the environment, depending on the shell you are using:

For sh or bash:
```
  cd ~/HCR_configuration/build
  source /set_build_env.sh
```  

For csh or tcsh:
```
  cd ~/HCR_configuration/build
  source /set_build_env.csh
```

This will set the following important environment variables:

```
 $LROSE_CORE_DIR: ~/git/lrose-core
 $LROSE_INSTALL_DIR: ~/lrose
 $RAP_MAKE_INC_DIR: include files used by the makefiles = ~/lrose/codebase/make_include
 $RAP_MAKE_BIN_DIR: scripts for the make = ~/lrose/codebase/make_bin
 $RAP_INC_DIR: the include install directory = ~/lrose/include
 $RAP_LIB_DIR: the library install directory = ~/lrose/lib
 $RAP_BIN_DIR: the binary install directory = ~/lrose/bin
 $RAP_SHARED_INC_DIR: the include install directory = ~/lrose/include
 $RAP_SHARED_LIB_DIR: the library install directory = ~/lrose/lib
 $RAP_SHARED_BIN_DIR: the binary install directory = ~/lrose/bin
```

All of these must be set for the build to work correctly.

After the HCR configuration is installed, these environment variables will be set
automatically via the ~/.cshrc file.

### Perform the build

To perform the build:
```
  cd ~/HCR_configuration/build
  source ./set_build_env.csh
  ./build_lrose.hcr
```

For the script details, see:

  See [build_lrose.hcr](../build/build_lrose.hcr)

