# determine mode of compilation
DEBUG  = YES
# choose solver
SOLVER = 
# name the compilation targets
TARGET = box interface

SPACEDIM=2
STRUCTURED=YES

APPFLAGS  = -DSPACEDIM=$(SPACEDIM)
ifeq ($(strip $(STRUCTURED)),YES)
	APPFLAGS += -DSTRUCTURED
endif

TESTMK=boxTest.mk interfaceTest.mk

# include configuration file
include $(INSILICOROOT)/config/convenience.mk

