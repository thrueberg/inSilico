# debugging flag
DEBUG=

# destination folder of binaries
BIN=$(INSILICOROOT)/tools/bin


# VTK headers and libraries
ifneq ($(origin VTK_PATH),undefined)
  VTK_INC = $(VTK_PATH)/include/vtk-5.8
  VTK_LIB_PATH = $(VTK_PATH)/lib/vtk-5.8
  VTK_LIB = -Wl,-rpath,$(VTK_LIB_PATH) -L$(VTK_LIB_PATH)
else
  VTK_PATH = /usr/include
  VTK_INC = $(VTK_PATH)/$(shell basename $(VTK_PATH)/vtk-5.*)
endif

INCLUDES = -I$(VTK_INC)
VTK_LIB  += -lvtkHybrid -lvtkIO -lvtkFiltering -lvtkCommon -lvtkGraphics

# libraries
LDLIBS  = $(VTK_LIB)

# ------------------------------------------------------------------------------
# standard targets, compilable with implicit rules
TARGET = vtk2vtu vtk2vts plotOverLine vtxDiff vtxComp

# ------------------------------------------------------------------------------
include $(INSILICOROOT)/config/convenience.mk

