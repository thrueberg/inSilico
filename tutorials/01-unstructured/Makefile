#######################################
# MAKEFILE for this directory         #
#######################################

include $(INSILICOROOT)/definitions.mk

# determine mode of compilation
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

BIN = $(HOME)/bin

# for cleanup
RM = -rm -f

TARGET = unstructured 


#######################################
# the rules
#
all: $(TARGET)

# clean: 
clean:
	$(RM) *.o $(TARGET) core

# This needs a more elegant solution
#
doc:
	$(BIN)/codeExtractor unstructured.cpp doc
	cd doc; asciidoc -a toc unstructured.txt

##doc: *.cpp
##	$(BIN)/codeExtractor $? doc

# implict rules:
# $(CC) $(LDFLAGS) n.o $(LOADLIBES) $(LDLIBS)
# $(CXX) -c $(CPPFLAGS) $(CXXFLAGS)