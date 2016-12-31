# Set environment variables for the make system
#
# for csh and tcsh
#

# HOST OS - most likely LINUX_64

setenv HOST_OS LINUX_64
uname -a | grep x86_64
if ($status == 1) then
    setenv HOST_OS LINUX
endif

# LROSE dirs

setenv LROSE_CORE_DIR ~/git/lrose-core
setenv LROSE_INSTALL_DIR ~/lrose

# make dirs

setenv RAP_MAKE_INC_DIR $LROSE_CORE_DIR/codebase/make_include
setenv RAP_MAKE_BIN_DIR $LROSE_CORE_DIR/codebase/make_bin

setenv RAP_INC_DIR $LROSE_INSTALL_DIR/include
setenv RAP_LIB_DIR $LROSE_INSTALL_DIR/lib
setenv RAP_BIN_DIR $LROSE_INSTALL_DIR/bin
setenv RAP_MAN_DIR $LROSE_INSTALL_DIR/man
setenv RAP_DOC_DIR $LROSE_INSTALL_DIR/doc

setenv RAP_SHARED_INC_DIR $LROSE_INSTALL_DIR/include
setenv RAP_SHARED_LIB_DIR $LROSE_INSTALL_DIR/lib
setenv RAP_SHARED_BIN_DIR $LROSE_INSTALL_DIR/bin
setenv RAP_SHARED_MAN_DIR $LROSE_INSTALL_DIR/man
setenv RAP_SHARED_DOC_DIR $LROSE_INSTALL_DIR/doc

setenv RAP_INST_LIB_DIR $LROSE_INSTALL_DIR/lib
setenv RAP_INST_BIN_DIR $LROSE_INSTALL_DIR/bin

# path

set path = ($RAP_BIN_DIR $path)

