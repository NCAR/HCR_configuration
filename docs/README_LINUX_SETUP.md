## Setting up LINUX for building LROSE

### gcc/g++ versions for LROSE build

LROSE expects support for the c++11 standard.

The gcc/g++ version should be 4.8.5 or later.

Distributions dated after June 2015 should work.

### Required LINUX and gcc/g++ versions for LROSE build

Most good, up-to date LINUX distributions should work.

Recommended distributions are:

  * RedHat
  * Centos
  * Fedora
  * Debian
  * Ubuntu

### Required LINUX packages for the LROSE build

For a full LROSE build under LINUX, you need the following packages:

```
  epel-release
  
  tcsh
  perl
  perl-Env

  ftp
  git
  emacs
  tkcvs

  gcc
  g++
  gfortran

  glibc-devel
  libX11-devel
  libXext-devel (if available)
  libpng-devel
  libtiff-devel
  jasper-devel
  zlib-devel
  libexpat-devel
  flex-devel
  fftw3-devel
  bzip2-devel
  jasper-devel
  qt4-devel

  gnuplot
  ImageMagick-devel
  ImageMagick-c++-devel

  xrdb
  Xvfb (virtual X server), specifically xorg-x11-server-Xvfb
  sshd (ssh logins)

  xorg-x11-fonts-misc
  xorg-x11-fonts-75dpi
  xorg-x11-fonts-100dpi
```

On Redhat-based hosts you can achieve this by running:

```
yum install -y epel-release \
tcsh perl perl-Env ftp git svn cvs tkcvs emacs tkcvs \
gcc gcc-c++ gcc-gfortran glibc-devel libX11-devel libXext-devel \
libpng-devel libtiff-devel jasper-devel zlib-devel expat-devel \
flex-devel fftw3-devel bzip2-devel jasper-devel qt4-devel xrdb \
Xvfb xorg-x11-fonts-misc xorg-x11-fonts-75dpi xorg-x11-fonts-100dpi \
gnuplot ImageMagick-devel ImageMagick-c++-devel
```

### Required LINUX packages for the CIDD build

For a full CIDD build under LINUX, you need the following packages:

```
  tcsh
  perl
  perl-Env

  ftp
  git

  gcc
  g++
  gfortran

  glibc-devel.i686
  libX11-devel.i686
  libXext-devel.i686 (if available)
  libtiff-devel.i686
  libpng-devel.i686
  libstdc++-devel.i686
  libtiff-devel.i686
  zlib-devel.i686
  expat-devel.i686
  flex-devel.i686
  fftw-devel.i686
  bzip2-devel.i686

  gnuplot
  ImageMagick-devel
  ImageMagick-c++-devel

  xrdb
  Xvfb (virtual X server), specifically xorg-x11-server-Xvfb
  sshd (ssh logins)

  xorg-x11-fonts-misc
  xorg-x11-fonts-75dpi
  xorg-x11-fonts-100dpi
```

On Redhat-based hosts you can achieve this by running:

```
yum install -y tcsh perl perl-Env ftp git svn cvs tkcvs emacs \
gcc gcc-c++ gcc-gfortran \
glibc-devel.i686 libX11-devel.i686 libXext-devel.i686 \
libtiff-devel.i686 libpng-devel.i686 \
libstdc++-devel.i686 libtiff-devel.i686 \
zlib-devel.i686 expat-devel.i686 flex-devel.i686 \
fftw-devel.i686 bzip2-devel.i686 xrdb Xvfb \
xorg-x11-fonts-misc xorg-x11-fonts-75dpi xorg-x11-fonts-100dpi \
gnuplot ImageMagick-devel ImageMagick-c++-devel
```

On Debian, you need to run the following:

```
  /usr/bin/dpkg --add-architecture i386
```

and use apt-get to install the following:

```
  libstdc++5:i386
  libstdc++6:i386
  libxml2:i386
  libgtk2.0-0:i386
  libgdk-pixbuf2.0-0:i386
```

### Running LROSE server-based applications

If you are making use of the data server applications in LROSE, you will need
to disable the firewall on your host, or open up the ports between 5300 and 5500.

To disable the firewall on a RedHat-based host:

```
  systemctl stop firewalld
  systemctl stop iptables
```

