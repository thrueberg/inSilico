# determine mode of compilation
DEBUG  = YES
# choose solver
SOLVER = SUPERLU
# name the compilation targets
TARGET = box interface hyperElastic

SPACEDIM=2
STRUCTURED=NO

APPFLAGS  = -DSPACEDIM=$(SPACEDIM)
ifeq ($(strip $(STRUCTURED)),YES)
	APPFLAGS += -DSTRUCTURED
endif

TESTMK=boxTest.mk interfaceTest.mk

# include configuration file
include $(INSILICOROOT)/config/convenience.mk

