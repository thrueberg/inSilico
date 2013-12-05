# ------------------------------------------------------------------------------
include $(INSILICOROOT)/definitions.mk
# ------------------------------------------------------------------------------
DEBUG=

# compile/link executables and flags
SHELL = /bin/bash
CXX   = $(CPLUSPLUS)
CC    = $(CPLUSPLUS)

# initialise includes
INCLUDES  = -I$(INSILICOROOT) -I$(BOOST_DIR)

# VTK headers and libraries
ifneq ($(origin VTK_PATH),undefined)
  VTK_INC = $(VTK_PATH)/include/vtk-5.8
  VTK_LIB_PATH = $(VTK_PATH)/lib/vtk-5.8
  VTK_LIB = -Wl,-rpath,$(VTK_LIB_PATH) -L$(VTK_LIB_PATH)
else
  VTK_PATH = /usr/include
  VTK_INC = $(VTK_PATH)/$(shell basename $(VTK_PATH)/vtk-5.*)
endif

INCLUDES += -I$(VTK_INC)
VTK_LIB  += -lvtkHybrid -lvtkIO -lvtkFiltering -lvtkCommon -lvtkGraphics

# compilation flags
CPPFLAGS = $(INCLUDES) -g3  #-03

# libraries
LIBS = $(VTK_LIB)


# set flags depending on mode
ifeq ($(DEBUG),YES) 
	CPPFLAGS += $(DEBDEFINES) $(DEBCPPFLAGS) 
	LDFLAGS  += $(DEBLDFLAGS)
else
	CPPFLAGS += $(RELDEFINES) $(RELCPPFLAGS) 
	LDFLAGS  += $(RELLDFLAGS)
endif

# libraries
LDLIBS  = $(LIBS)

# bin folder
ifneq (,$(wildcard $(HOME)/bin))
	BIN = $(HOME)/bin
else
	BIN = /dev/null
endif

# ------------------------------------------------------------------------------
# standard targets, compilable with implicit rules
TARGETS = vtk2vtu

# ------------------------------------------------------------------------------
all: $(TARGETS) 

clean:
	rm -f  *.o $(TARGETS) core

install: $(TARGETS)
	@echo make install: Copying $(TARGETS) to executable folder $(BIN)
	@cp $(TARGETS) $(BIN)

