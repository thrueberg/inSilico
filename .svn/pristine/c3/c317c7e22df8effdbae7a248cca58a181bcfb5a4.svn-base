#######################################
# MAKEFILE for this directory         #
#######################################

include $(INSILICOROOT)/definitions.mk

# determine modus of compilation
DEBUG=YES
SOLVER = 

# compile/link executables and flags
SHELL = /bin/bash
CXX   = $(CPLUSPLUS)
CC    = $(CPLUSPLUS)

# includes
INCLUDES  = -I$(INSILICOROOT) 
INCLUDES += -I$(BOOST_DIR)
INCLUDES += -isystem$(EIGEN3_DIR)  #note the -isystem for surpression of warnings

# libraries
LDLIBS  = 

# common flags
CPPFLAGS = $(INCLUDES) 

# set flags depending on mode
ifeq ($(strip $(DEBUG)),YES) 
	CPPFLAGS += $(DEBDEFINES) $(DEBCPPFLAGS) 
	LDFLAGS  += $(DEBLDFLAGS)
else
	CPPFLAGS += $(RELDEFINES) $(RELCPPFLAGS) 
	LDFLAGS  += $(RELLDFLAGS)
endif

# bin folder
ifneq (,$(wildcard $(HOME)/bin))
	BIN = $(HOME)/bin
else
	BIN = /dev/null
endif


# for cleanup
RM = -rm -f

TARGETS = smfReverse

#######################################
# the rules
#
all: $(TARGETS)

# clean: 
clean:
	$(RM) *.o $(TARGETS) core

install: $(TARGETS)
	@echo make install: Copying $(TARGETS) to executable folder $(BIN)
	@cp $(TARGETS) $(BIN)


# implict rules:
# $(CC) $(LDFLAGS) n.o $(LOADLIBES) $(LDLIBS)
# $(CXX) -c $(CPPFLAGS) $(CXXFLAGS)