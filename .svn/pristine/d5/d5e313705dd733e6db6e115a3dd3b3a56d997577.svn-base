################################################################################
# Include machine dependent settings file
include $(INSILICOROOT)/config/definitions.mk

# default modes
DEBUG    ?= YES
SOLVER   ?= EIGEN
APPFLAGS ?=
NTHREADS ?= 1
VTK      ?= NO

# compile/link executables and flags
SHELL = /bin/bash
CXX   = $(CPLUSPLUS)
CC    = $(CPLUSPLUS)

# includes
INCLUDES ?=
INCLUDES += $(SYSINCLUDES)

# libraries
LDLIBS  ?= 

# common flags
CPPFLAGS = $(APPFLAGS) $(INCLUDES) -DNTHREADS=$(NTHREADS)

################################################################################
# set flags depending on mode (DEBUG or RELEASE)
ifeq ($(strip $(DEBUG)),YES) 
	CPPFLAGS += $(DEBDEFINES) $(DEBCPPFLAGS) 
	LDFLAGS  += $(DEBLDFLAGS)
else
	CPPFLAGS += $(RELDEFINES) $(RELCPPFLAGS) 
	LDFLAGS  += $(RELLDFLAGS)
endif

################################################################################
# set flags depending on solver choice

# SUPERLU
ifeq ($(findstring SUPERLU,$(SOLVER)),SUPERLU)
	CPPFLAGS += -DLOAD_SUPERLU
	INCLUDES += -I$(SUPERLU_INC)
	LDLIBS   += -L $(SUPERLU_LIB_PATH) $(SUPERLU_LIB)
endif

# PARDISO
ifeq ($(findstring PARDISO,$(SOLVER)),PARDISO)
	CPPFLAGS += -DLOAD_PARDISO
  ifneq ($(INSILICOCOMPENV),INTEL)
	REQUIREINTEL=Pardiso
  endif
endif

# UMFPACK
ifeq ($(findstring UMFPACK,$(SOLVER)),UMFPACK)
	CPPFLAGS += -DLOAD_UMFPACK
	INCLUDES += -I$(UMFPACK_INC)
	LDLIBS   += -L $(UMFPACK_LIB_PATH) $(UMFPACK_LIB)
endif

# VTK
ifeq ($(strip $(VTK)),YES)
	INCLUDES += -I$(VTK_INC)
	LDLIBS   +=   $(VTK_LIB)
endif


# error message
ifdef REQUIREINTEL
$(error  (EE) Intel compilation environment required by $(REQUIREINTEL))
endif

# for cleanup
RM = -rm -f

################################################################################
# test scripts (if not set)
TESTMK ?= appTest.mk

# if test make file is not available, write a message
ifeq ($(wildcard $(TESTMK)),) 
	APPTEST=@echo No test makefile $(TESTMK) available
else
	APPTEST=($(foreach BLA,$(TESTMK),\
	(make -j1 -f $(BLA) all && make -j1 -f $(BLA) triumph) &&) true) || \
	(make -j1 -f $(TESTMK) defeat && false) 
endif
# make -j1 -f $(TESTMK) all && make -j1 -f $(TESTMK) triumph)|| \

################################################################################
INSTALL=
ifneq ($(BIN),)
	INSTALL=@cp $(TARGET) $(BIN)
endif

################################################################################
# the rules
#
all: $(TARGET)

# clean: 
clean:
	$(RM) *.o $(TARGET) core

# copy binaries to bin folder
install: $(TARGET) 
	$(INSTALL)

# call application test funtionality
appTest:
	$(APPTEST)


# implict rules:
# $(CC) $(LDFLAGS) n.o $(LOADLIBES) $(LDLIBS)
# $(CXX) -c $(CPPFLAGS) $(CXXFLAGS)
